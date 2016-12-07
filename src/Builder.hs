-- | Buffer builder to assemble Bgzf blocks.  The plan is to serialize
-- stuff (BAM and BCF) into a buffer, then Bgzf chunks from the buffer.
-- We use a large buffer, and we always make sure there is plenty of
-- space in it (to avoid redundant checks).  Whenever a block is ready
-- to be compressed, we stick it into a MVar.  When we run out of space,
-- we simply use a new buffer.  Multiple threads grab pieces from the
-- MVar, compress them, pass them downstream through another MVar.  A
-- final thread restores the order and writes the blocks.

module BgzfBuilder where

-- import Bio.Iteratee hiding ( NullPoint )
-- import Bio.Iteratee.Bgzf
import Bio.Prelude
import Data.NullPoint ( NullPoint(..) )
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable

import qualified Data.ByteString            as B
import qualified Data.ByteString.Unsafe     as B
import qualified Data.ByteString.Builder    as B ( Builder, toLazyByteString )
import qualified Data.ByteString.Lazy       as B ( foldrChunks )

-- | The 'MutableByteArray' is garbage collected, so we don't get leaks.
-- Once it has grown to a practical size (and the initial 128k should be
-- very practical), we don't get fragmentation either.  We also avoid
-- copies for the most part, since no intermediate 'ByteString's, either
-- lazy or strict have to be allocated.
data BB = BB { buffer :: {-# UNPACK #-} !(ForeignPtr Word8)
             , size   :: {-# UNPACK #-} !Int
             , len    :: {-# UNPACK #-} !Int
             , mark   :: {-# UNPACK #-} !Int
             , mark2  :: {-# UNPACK #-} !Int }

-- This still seems to have considerable overhead.  Don't know if this
-- can be improved by effectively inlining IO and turning the BB into an
-- unboxed tuple.  XXX
newtype Push = Push (BB -> IO BB)

instance Monoid Push where
    {-# INLINE mempty #-}
    mempty                  = Push return
    {-# INLINE mappend #-}
    Push a `mappend` Push b = Push (a >=> b)

instance NullPoint Push where
    empty = Push return


-- | Creates a buffer with a given initial capacity.
newBuffer :: Int -> IO BB
newBuffer sz = mallocForeignPtrBytes sz >>= \ar -> return $ BB ar sz 0 0 0

-- | Ensures a given free space in the buffer by doubling its capacity
-- if necessary.
{-# INLINE ensureBuffer #-}
ensureBuffer :: Int -> Push
ensureBuffer n = Push $ \b ->
    if len b + n < size b
    then return b
    else expandBuffer b

expandBuffer :: BB -> IO BB
expandBuffer b = do arr1 <- mallocForeignPtrBytes (size b + size b)
                    withForeignPtr arr1 $ \d ->
                        withForeignPtr (buffer b) $ \s ->
                             copyBytes d s (len b)
                    return $ BB { buffer = arr1
                                , size   = size b + size b
                                , len    = len b
                                , mark   = mark b
                                , mark2  = mark2 b }

{-# INLINE unsafePushByte #-}
unsafePushByte :: Word8 -> Push
unsafePushByte w = Push $ \b -> do
    withForeignPtr (buffer b) $ \p ->
        pokeByteOff p (len b) w
    return $ b { len = len b + 1 }

{-# INLINE pushByte #-}
pushByte :: Word8 -> Push
pushByte b = ensureBuffer 1 <> unsafePushByte b

{-# INLINE unsafePushWord32 #-}
unsafePushWord32 :: Word32 -> Push
unsafePushWord32 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 16)
                  <> unsafePushByte (fromIntegral $ w `shiftR` 24)

{-# INLINE unsafePushWord16 #-}
unsafePushWord16 :: Word16 -> Push
unsafePushWord16 w = unsafePushByte (fromIntegral $ w `shiftR`  0)
                  <> unsafePushByte (fromIntegral $ w `shiftR`  8)

{-# INLINE pushWord32 #-}
pushWord32 :: Word32 -> Push
pushWord32 w = ensureBuffer 4 <> unsafePushWord32 w

{-# INLINE pushWord16 #-}
pushWord16 :: Word16 -> Push
pushWord16 w = ensureBuffer 2 <> unsafePushWord16 w

{-# INLINE unsafePushByteString #-}
unsafePushByteString :: B.ByteString -> Push
unsafePushByteString bs = Push $ \b ->
    B.unsafeUseAsCStringLen bs $ \(p,ln) ->
        withForeignPtr (buffer b)  $ \adr ->
            b { len = len b + ln } <$
                copyBytes (adr `plusPtr` len b) p ln

{-# INLINE pushByteString #-}
pushByteString :: B.ByteString -> Push
pushByteString bs = ensureBuffer (B.length bs) <> unsafePushByteString bs

{-# INLINE unsafePushFloat #-}
unsafePushFloat :: Float -> Push
unsafePushFloat f =
    unsafePushWord32 $ unsafeDupablePerformIO $
    alloca $ \b -> poke (castPtr b) f >> peek b

{-# INLINE pushFloat #-}
pushFloat :: Float -> Push
pushFloat f = ensureBuffer 4 <> unsafePushFloat f

{-# INLINE pushBuilder #-}
pushBuilder :: B.Builder -> Push
pushBuilder = B.foldrChunks ((<>) . pushByteString) mempty . B.toLazyByteString

-- | Sets a mark.  This can later be filled in with a record length
-- (used to create BAM records).
{-# INLINE unsafeSetMark #-}
unsafeSetMark :: Push
unsafeSetMark = Push $ \b -> return $ b { len = len b + 4, mark = len b }

{-# INLINE setMark #-}
setMark :: Push
setMark = ensureBuffer 4 <> unsafeSetMark

-- | Ends a record by filling the length into the field that was
-- previously marked.  Terrible things will happen if this wasn't
-- preceded by a corresponding 'setMark'.
{-# INLINE endRecord #-}
endRecord :: Push
endRecord = Push $ \b -> withForeignPtr (buffer b) $ \p -> do
    let !l = len b - mark b - 4
    pokeByteOff p (mark b + 0) (fromIntegral $ shiftR l  0 :: Word8)
    pokeByteOff p (mark b + 1) (fromIntegral $ shiftR l  8 :: Word8)
    pokeByteOff p (mark b + 2) (fromIntegral $ shiftR l 16 :: Word8)
    pokeByteOff p (mark b + 3) (fromIntegral $ shiftR l 24 :: Word8)
    return b

-- | Ends the first part of a record.  The length is filled in *before*
-- the mark, which is specifically done to support the *two* length
-- fields in BCF.  It also remembers the current position.  Horrible
-- things happen if this isn't preceeded by *two* succesive invocations
-- of 'setMark'.
{-# INLINE endRecordPart1 #-}
endRecordPart1 :: Push
endRecordPart1 = Push $ \b -> withForeignPtr (buffer b) $ \p -> do
    let !l = len b - mark b - 4
    pokeByteOff p (mark b - 4) (fromIntegral $ shiftR l  0 :: Word8)
    pokeByteOff p (mark b - 3) (fromIntegral $ shiftR l  8 :: Word8)
    pokeByteOff p (mark b - 2) (fromIntegral $ shiftR l 16 :: Word8)
    pokeByteOff p (mark b - 1) (fromIntegral $ shiftR l 24 :: Word8)
    return $ b { mark2 = len b }

-- | Ends the second part of a record.  The length is filled in at the
-- mark, but computed from the sencond mark only.  This is specifically
-- done to support the *two* length fields in BCF.  Horrible things
-- happen if this isn't preceeded by *two* succesive invocations of
-- 'setMark' and one of 'endRecordPart1'.
{-# INLINE endRecordPart2 #-}
endRecordPart2 :: Push
endRecordPart2 = Push $ \b -> withForeignPtr (buffer b) $ \p -> do
    let !l = len b - mark2 b
    pokeByteOff p (mark b + 0) (fromIntegral $ shiftR l  0 :: Word8)
    pokeByteOff p (mark b + 1) (fromIntegral $ shiftR l  8 :: Word8)
    pokeByteOff p (mark b + 2) (fromIntegral $ shiftR l 16 :: Word8)
    pokeByteOff p (mark b + 3) (fromIntegral $ shiftR l 24 :: Word8)
    return b


{-# INLINE encodeBgzfWith #-}
encodeBgzfWith :: MonadIO m => Int -> Enumeratee Push B.ByteString m b
encodeBgzfWith lv o = newBuffer 128000 `ioBind` \bb -> eneeCheckIfDone (liftI . step bb) o
  where
    step bb k (EOF  mx) = finalFlush bb k mx
    step bb k (Chunk (Push p)) = p bb `ioBind` \bb' -> tryFlush bb' 0 k

    tryFlush bb off k
        | len bb - off < maxBlockSize
            = withForeignPtr (buffer bb)
                    (\p -> moveBytes p (p `plusPtr` off) (len bb - off))
              `ioBind_` liftI (step (bb { len = len bb - off
                                        , mark = mark bb - off `max` 0 }) k)
        | otherwise
            = withForeignPtr (buffer bb)
                    (\adr -> compressChunk lv (adr `plusPtr` off) (fromIntegral maxBlockSize))
              `ioBind` eneeCheckIfDone (tryFlush bb (off+maxBlockSize)) . k . Chunk

    finalFlush bb k mx
        | len bb < maxBlockSize
            = withForeignPtr (buffer bb)
                    (\adr -> compressChunk lv (castPtr adr) (fromIntegral $ len bb))
              `ioBind` eneeCheckIfDone (finalFlush2 mx) . k . Chunk

        | otherwise
            = error "WTF?!  This wasn't supposed to happen."

    finalFlush2 mx k = idone (k $ Chunk bgzfEofMarker) (EOF mx)



