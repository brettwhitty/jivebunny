-- | Handling of BGZF files.  This time, there won't be a dependency on
-- Iteratee.

module Bgzf (
    Block(..), decompressBgzfBlocks', decompressBgzfBlocks,
    decompressBgzf, decompressPlain,
    maxBlockSize, bgzfEofMarker, liftBlock, getOffset,
    BgzfChunk(..), isBgzf, isGzip, parMapChunksIO,
    compressBgzf, compressBgzfLv, compressBgzf', CompressParams(..),
    compressChunk
                     ) where

import Bio.Iteratee
import Bio.Prelude
import Control.Concurrent.Async             ( async, wait )
import Foreign.Marshal.Alloc                ( mallocBytes, free, allocaBytes )
import Foreign.Storable                     ( peekByteOff, pokeByteOff )
import Foreign.C.String                     ( withCAString )
import Foreign.C.Types                      ( CInt(..), CChar(..), CUInt(..), CULong(..) )
import Foreign.Ptr                          ( nullPtr, castPtr, Ptr, plusPtr, minusPtr )

import qualified Data.ByteString            as S
import qualified Data.ByteString.Unsafe     as S
import qualified Data.Iteratee.ListLike     as I

#include <zlib.h>

-- | One BGZF block: virtual offset and contents.  Could also be a block
-- of an uncompressed file, if we want to support indexing of
-- uncompressed BAM or some silliness like that.
data Block = Block { block_offset   :: {-# UNPACK #-} !FileOffset
                   , block_contents :: {-# UNPACK #-} !Bytes }

instance NullPoint Block where empty = mempty
instance Nullable Block where nullC (Block _ s) = S.null s

instance Monoid Block where
    mempty = Block 0 S.empty
    mappend (Block x s) (Block _ t) = Block x (s `S.append` t)
    mconcat [] = Block 0 S.empty
    mconcat bs@(Block x _:_) = Block x $ S.concat [s|Block _ s <- bs]

-- | "Decompresses" a plain file.  What's actually happening is that the
-- offset in the input stream is tracked and added to the @Bytes@s
-- giving @Block@s.  This results in the same interface as decompressing
-- actual Bgzf.
decompressPlain :: MonadIO m => Enumeratee Bytes Block m a
decompressPlain = eneeCheckIfDone (liftI . step 0)
  where
    step !o it (Chunk s) = eneeCheckIfDone (liftI . step (o + fromIntegral (S.length s))) . it $ Chunk (Block o s)
    step  _ it (EOF  mx) = idone (liftI it) (EOF mx)

-- | Decompress a BGZF stream into a stream of 'Bytes's.
decompressBgzf :: MonadIO m => Enumeratee Bytes Bytes m a
decompressBgzf = decompressBgzfBlocks ><> mapChunks block_contents

decompressBgzfBlocks :: MonadIO m => Enumeratee Bytes Block m a
decompressBgzfBlocks out =  do
    np <- liftIO $ getNumCapabilities
    decompressBgzfBlocks' np out

-- | Decompress a BGZF stream into a stream of 'Block's, 'np' fold parallel.
decompressBgzfBlocks' :: MonadIO m => Int -> Enumeratee Bytes Block m a
decompressBgzfBlocks' np = eneeCheckIfDonePass (go 0 emptyQ)
  where
    -- check if the queue is full
    go !off !qq k (Just e) = handleSeek off qq k e
    go !off !qq k Nothing = case popQ qq of
        Just (a, qq') | lengthQ qq == np -> liftIO (wait a) >>= eneeCheckIfDonePass (go off qq') . k . Chunk
        _                                -> liftI $ go' off qq k

    -- we have room for input, so try and get a compressed block
    go' !_   !qq k (EOF  mx) = goE mx qq k Nothing
    go' !off !qq k (Chunk c)
        | S.null  c = liftI $ go' off qq k
        | otherwise = joinIM $ enumPure1Chunk c $ do
                                  (off', op) <- get_bgzf_block off
                                  a <- liftIO (async op)
                                  go off' (pushQ a qq) k Nothing

    -- input ended, empty the queue
    goE  _ !qq k (Just e) = handleSeek 0 qq k e
    goE mx !qq k Nothing = case popQ qq of
        Nothing      -> idone (liftI k) (EOF mx)
        Just (a,qq') -> liftIO (wait a) >>= eneeCheckIfDonePass (goE mx qq') . k . Chunk

    handleSeek !off !qq k e = case fromException e of
        Nothing                -> throwRecoverableErr e $ go' off qq k
        Just (SeekException o) -> do
            cancelAll qq
            seek $ o `shiftR` 16
            eneeCheckIfDonePass (go (o `shiftR` 16) emptyQ) $ do
                block'drop . fromIntegral $ o .&. 0xffff
                k (EOF Nothing)
                -- I think, 'seek' swallows one 'Stream' value on
                -- purpose, so we have to give it a dummy one.

    block'drop sz = liftI $ \s -> case s of
        EOF _ -> throwErr $ setEOF s
        Chunk (Block p c)
            | S.length c < sz -> block'drop (sz - S.length c)
            | otherwise       -> let b' = Block (p + fromIntegral sz) (S.drop sz c)
                                 in idone () (Chunk b')

get_bgzf_block :: MonadIO m => FileOffset -> Iteratee Bytes m (FileOffset, IO Block)
get_bgzf_block off = do !(csize,xlen) <- get_bgzf_header
                        !comp  <- get_block . fromIntegral $ csize - xlen - 19
                        !crc   <- endianRead4 LSB
                        !isize <- endianRead4 LSB

                        let !off' = off + fromIntegral csize + 1
                            op    = decompress1 (off `shiftL` 16) comp crc (fromIntegral isize)
                        return (off',op)
  where
    -- Get a block of a prescribed size.  Comes back as a list of chunks.
    get_block sz = liftI $ \s -> case s of
        EOF _ -> throwErr $ setEOF s
        Chunk c | S.length c < sz -> (:) c `liftM` get_block (sz - S.length c)
                | otherwise       -> idone [S.take sz c] (Chunk (S.drop sz c))


-- | Decodes a BGZF block header and returns the block size if
-- successful.
get_bgzf_header :: Monad m => Iteratee Bytes m (Word16, Word16)
get_bgzf_header = do n <- I.heads "\31\139"
                     _cm <- I.head
                     flg <- I.head
                     if flg `testBit` 2 then do
                         I.drop 6
                         xlen <- endianRead2 LSB
                         it <- I.take (fromIntegral xlen) get_bsize >>= lift . tryRun
                         case it of Left e -> throwErr e
                                    Right s | n == 2 -> return (s,xlen)
                                    _ -> throwErr $ iterStrExc "No BGZF"
                      else throwErr $ iterStrExc "No BGZF"
  where
    get_bsize = do i1 <- I.head
                   i2 <- I.head
                   len <- endianRead2 LSB
                   if i1 == 66 && i2 == 67 && len == 2
                      then endianRead2 LSB
                      else I.drop (fromIntegral len) >> get_bsize

-- | Tests whether a stream is in BGZF format.  Does not consume any
-- input.
isBgzf :: Monad m => Iteratee Bytes m Bool
isBgzf = liftM isRight $ checkErr $ iLookAhead $ get_bgzf_header

-- | Tests whether a stream is in GZip format.  Also returns @True@ on a
-- Bgzf stream, which is technically a special case of GZip.
isGzip :: Monad m => Iteratee Bytes m Bool
isGzip = liftM (either (const False) id) $ checkErr $ iLookAhead $ test
  where
    test = do n <- I.heads "\31\139"
              I.drop 24
              b <- I.isFinished
              return $ not b && n == 2

-- ------------------------------------------------------------------------- Output

-- | Maximum block size for Bgzf: 64k with some room for headers and
-- uncompressible stuff
maxBlockSize :: Int
maxBlockSize = 65450


-- | The EOF marker for BGZF files.
-- This is just an empty string compressed as BGZF.  Appended to BAM
-- files to indicate their end.
bgzfEofMarker :: Bytes
bgzfEofMarker = "\x1f\x8b\x8\x4\0\0\0\0\0\xff\x6\0\x42\x43\x2\0\x1b\0\x3\0\0\0\0\0\0\0\0\0"

-- | Decompress a collection of strings into a single BGZF block.
--
-- Ideally, we receive one decode chunk from a BGZF file, decompress it,
-- and return it, in the process attaching the virtual address.  But we
-- might actually get more than one chunk, depending on the internals of
-- the @Iteratee@s used.  If so, we concatenate them; the first gets to
-- assign the address.
--
-- Now allocate space for uncompressed data, decompress the chunks we
-- got, compute crc for each and check it, finally convert to 'Bytes'
-- and emit.
--
-- We could probably get away with @unsafePerformIO@'ing everything in
-- here, but then again, we only do this when we're writing output
-- anyway.  Hence, run in IO.


decompress1 :: FileOffset -> [Bytes] -> Word32 -> Int -> IO Block
decompress1 off ss crc usize =
    allocaBytes (#{const sizeof(z_stream)}) $ \stream -> do
    buf <- mallocBytes usize

    #{poke z_stream, msg}       stream nullPtr
    #{poke z_stream, zalloc}    stream nullPtr
    #{poke z_stream, zfree}     stream nullPtr
    #{poke z_stream, opaque}    stream nullPtr
    #{poke z_stream, next_in}   stream nullPtr
    #{poke z_stream, next_out}  stream buf
    #{poke z_stream, avail_in}  stream (0 :: CUInt)
    #{poke z_stream, avail_out} stream (fromIntegral usize :: CUInt)

    z_check "inflateInit2" =<< c_inflateInit2 stream (-15)

    -- loop over the fragments, forward order
    forM_ ss $ \s -> case fromIntegral $ S.length s of
            l | l > 0 -> S.unsafeUseAsCString s $ \p -> do
                #{poke z_stream, next_in} stream p
                #{poke z_stream, avail_in} stream (l :: CUInt)
                z_check "inflate" =<< c_inflate stream #{const Z_NO_FLUSH}
            _ -> return ()

    z_check "inflate" =<< c_inflate stream #{const Z_FINISH}
    z_check "inflateEnd" =<< c_inflateEnd stream

    pe <- #{peek z_stream, next_out} stream
    when (pe `minusPtr` buf /= usize) $ error "size mismatch after deflate()"

    crc0 <- c_crc32 0 nullPtr 0
    crc' <- c_crc32 crc0 buf (fromIntegral usize)
    when (fromIntegral crc /= crc') $ error "CRC error after deflate()"

    Block off `liftM` S.unsafePackCStringFinalizer (castPtr buf) usize (free buf)


-- | Compress a collection of strings into a single BGZF block.
--
-- Okay, performance was lacking... let's do it again, in a more direct
-- style.  We build our block manually.  First check if the compressed
-- data is going to fit---if not, that's a bug.  Then alloc a buffer,
-- fill with a dummy header, alloc a ZStream, compress the pieces we
-- were handed one at a time.  Calculate CRC32, finalize header,
-- construct a byte string, return it.
--
-- We could probably get away with @unsafePerformIO@'ing everything in
-- here, but then again, we only do this when we're writing output
-- anyway.  Hence, run in IO.

compress1 :: Int -> [Bytes] -> IO Bytes
compress1 _lv [] = return bgzfEofMarker
compress1 lv ss0 =
    allocaBytes (#{const sizeof(z_stream)}) $ \stream -> do

    let input_length = sum (map S.length ss0)
    when (input_length > maxBlockSize) $ error "Trying to create too big a BGZF block; this is a bug."
    buf <- mallocBytes 65536

    -- steal header from the EOF marker (length is wrong for now)
    S.unsafeUseAsCString bgzfEofMarker $ \eof ->
        forM_ [0,4..16] $ \o -> do x <- peekByteOff eof o
                                   pokeByteOff buf o (x::Word32)

    #{poke z_stream, msg}       stream nullPtr
    #{poke z_stream, zalloc}    stream nullPtr
    #{poke z_stream, zfree}     stream nullPtr
    #{poke z_stream, opaque}    stream nullPtr
    #{poke z_stream, next_in}   stream nullPtr
    #{poke z_stream, next_out}  stream (buf `plusPtr` 18)
    #{poke z_stream, avail_in}  stream (0 :: CUInt)
    #{poke z_stream, avail_out} stream (65536-18-8 :: CUInt)

    z_check "deflateInit2" =<< c_deflateInit2 stream (fromIntegral lv) #{const Z_DEFLATED}
                                              (-15) 8 #{const Z_DEFAULT_STRATEGY}

    -- loop over the fragments.  In reverse order!
    let go (s:ss) = do
            crc <- go ss
            S.unsafeUseAsCString s $ \p ->
              case fromIntegral $ S.length s of
                l | l > 0 -> do
                    #{poke z_stream, next_in} stream p
                    #{poke z_stream, avail_in} stream (l :: CUInt)
                    z_check "deflate" =<< c_deflate stream #{const Z_NO_FLUSH}
                    c_crc32 crc p l
                _ -> return crc
        go [] = c_crc32 0 nullPtr 0
    crc <- go ss0

    z_check "deflate" =<< c_deflate stream #{const Z_FINISH}
    z_check "deflateEnd" =<< c_deflateEnd stream

    compressed_length <- (+) (18+8) `fmap` #{peek z_stream, total_out} stream
    when (compressed_length > 65536) $ error "produced too big a block"

    -- set length in header
    pokeByteOff buf 16 (fromIntegral $ (compressed_length-1) .&. 0xff :: Word8)
    pokeByteOff buf 17 (fromIntegral $ (compressed_length-1) `shiftR` 8 :: Word8)

    pokeByteOff buf (compressed_length-8) (fromIntegral crc :: Word32)
    pokeByteOff buf (compressed_length-4) (fromIntegral input_length :: Word32)

    S.unsafePackCStringFinalizer buf compressed_length (free buf)


data ZStream

{-# INLINE z_check #-}
z_check :: String -> CInt -> IO ()
z_check msg c = when (c /= #{const Z_OK} && c /= #{const Z_STREAM_END}) $
                   error $ msg ++ " failed: " ++ show c


c_deflateInit2 :: Ptr ZStream -> CInt -> CInt -> CInt -> CInt -> CInt -> IO CInt
c_deflateInit2 z a b c d e = withCAString #{const_str ZLIB_VERSION} $ \versionStr ->
    c_deflateInit2_ z a b c d e versionStr (#{const sizeof(z_stream)} :: CInt)

foreign import ccall unsafe "zlib.h deflateInit2_" c_deflateInit2_ ::
    Ptr ZStream -> CInt -> CInt -> CInt -> CInt -> CInt
                -> Ptr CChar -> CInt -> IO CInt

c_inflateInit2 :: Ptr ZStream -> CInt -> IO CInt
c_inflateInit2 z a = withCAString #{const_str ZLIB_VERSION} $ \versionStr ->
    c_inflateInit2_ z a versionStr (#{const sizeof(z_stream)} :: CInt)

foreign import ccall unsafe "zlib.h inflateInit2_" c_inflateInit2_ ::
    Ptr ZStream -> CInt -> Ptr CChar -> CInt -> IO CInt

foreign import ccall unsafe "zlib.h deflate" c_deflate ::
    Ptr ZStream -> CInt -> IO CInt

foreign import ccall unsafe "zlib.h inflate" c_inflate ::
    Ptr ZStream -> CInt -> IO CInt

foreign import ccall unsafe "zlib.h deflateEnd" c_deflateEnd ::
    Ptr ZStream -> IO CInt

foreign import ccall unsafe "zlib.h inflateEnd" c_inflateEnd ::
    Ptr ZStream -> IO CInt

foreign import ccall unsafe "zlib.h crc32" c_crc32 ::
    CULong -> Ptr CChar -> CUInt -> IO CULong

-- ------------------------------------------------------------------------------------------------- utils

-- | Get the current virtual offset.  The virtual address in a BGZF
-- stream contains the offset of the current block in the upper 48 bits
-- and the current offset into that block in the lower 16 bits.  This
-- scheme is compatible with the way BAM files are indexed.
getOffset :: Iteratee Block m FileOffset
getOffset = liftI step
  where
    step s@(EOF _) = icont step (Just (setEOF s))
    step s@(Chunk (Block o _)) = idone o s

-- | Runs an @Iteratee@ for @Bytes@s when decompressing BGZF.  Adds
-- internal bookkeeping.
liftBlock :: Monad m => Iteratee Bytes m a -> Iteratee Block m a
liftBlock = liftI . step
  where
    step it (EOF ex) = joinI $ lift $ enumChunk (EOF ex) it

    step it (Chunk (Block !l !s)) = Iteratee $ \od oc ->
            enumPure1Chunk s it >>= \it' -> runIter it' (onDone od) (oc . step . liftI)
      where
        !sl = S.length s
        onDone od hdr (Chunk !rest) = od hdr . Chunk $! Block (l + fromIntegral (sl - S.length rest)) rest
        onDone od hdr (EOF      ex) = od hdr (EOF ex)


-- | Compresses a stream of @Bytes@s into a stream of BGZF blocks,
-- in parallel

-- We accumulate an uncompressed block as long as adding a new chunk to
-- it doesn't exceed the max. block size.  If we receive an empty chunk
-- (used as a flush signal), or if we would exceed the block size, we
-- write out a block.  Then we continue writing until we're below block
-- size.  On EOF, we flush and write the end marker.

compressBgzf' :: MonadIO m => CompressParams -> Enumeratee BgzfChunk Bytes m a
compressBgzf' (CompressParams lv np) = bgzfBlocks ><> parMapChunksIO np (compress1 lv)

data BgzfChunk = SpecialChunk  !Bytes BgzfChunk
               | RecordChunk   !Bytes BgzfChunk
               | LeftoverChunk !Bytes BgzfChunk
               | NoChunk

instance NullPoint BgzfChunk where empty = NoChunk
instance Nullable BgzfChunk where
    nullC NoChunk = True
    nullC (SpecialChunk  s c) = S.null s && nullC c
    nullC (RecordChunk   s c) = S.null s && nullC c
    nullC (LeftoverChunk s c) = S.null s && nullC c

-- | Breaks a stream into chunks suitable to be compressed individually.
-- Each chunk on output is represented as a list of 'Bytes',
-- each list must be reversed and concatenated to be compressed.
-- ('compress1' does that.)

bgzfBlocks :: Monad m => Enumeratee BgzfChunk [Bytes] m a
bgzfBlocks = eneeCheckIfDone (liftI . to_blocks 0 [])
  where
    -- terminate by sending the last block and then an empty block,
    -- which becomes the EOF marker
    to_blocks _alen acc k (EOF mx) =
        lift (enumPure1Chunk [S.empty] (k $ Chunk acc)) >>= flip idone (EOF mx)

    -- \'Empty list\', in a sense.
    to_blocks  alen acc k (Chunk NoChunk) = liftI $ to_blocks alen acc k

    to_blocks  alen acc k (Chunk (SpecialChunk c cs))  -- special chunk, encode then flush
        -- If it fits, flush.
        | alen + S.length c < maxBlockSize  = eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk cs)) . k $ Chunk (c:acc)
        -- If nothing is pending, flush the biggest thing that does fit.
        | null acc                       = let (l,r) = S.splitAt maxBlockSize c
                                           in eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (SpecialChunk r cs))) . k $ Chunk [l]
        -- Otherwise, flush what's pending and think again.
        | otherwise                         = eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (SpecialChunk c cs))) . k $ Chunk acc

    to_blocks  alen acc k (Chunk (RecordChunk c cs))
        -- if it fits, we accumulate,  (needs to consider the length prefix!)
        | alen + S.length c + 4 < maxBlockSize  = to_blocks (alen + S.length c + 4) (c:encLength c:acc) k (Chunk cs)
        -- else if nothing's pending, we break the chunk,  (needs to consider the length prefix!)
        | null acc                       = let (l,r) = S.splitAt (maxBlockSize-4) c
                                           in eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (LeftoverChunk r cs))) . k $
                                                    Chunk [l, encLength l]
        -- else we flush the accumulator and think again.
        | otherwise                         = eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (RecordChunk c cs))) . k $ Chunk acc
      where
        encLength s = let !l = S.length s in S.pack [ fromIntegral (l `shiftR`  0 .&. 0xff)
                                                    , fromIntegral (l `shiftR`  8 .&. 0xff)
                                                    , fromIntegral (l `shiftR` 16 .&. 0xff)
                                                    , fromIntegral (l `shiftR` 24 .&. 0xff) ]

    to_blocks  alen acc k (Chunk (LeftoverChunk c cs))
        -- if it fits, we accumulate,
        | alen + S.length c < maxBlockSize  = to_blocks (alen + S.length c) (c:acc) k (Chunk cs)
        -- else if nothing's pending, we break the chunk,
        | null acc                       = let (l,r) = S.splitAt maxBlockSize c
                                           in eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (LeftoverChunk r cs))) . k $ Chunk [l]
        -- else we flush the accumulator and think again.
        | otherwise                         = eneeCheckIfDone (\k' -> to_blocks 0 [] k' (Chunk (LeftoverChunk c cs))) . k $ Chunk acc

-- | Like 'compressBgzf'', with sensible defaults.
compressBgzf :: MonadIO m => Enumeratee BgzfChunk Bytes m a
compressBgzf = compressBgzfLv 6

compressBgzfLv :: MonadIO m => Int -> Enumeratee BgzfChunk Bytes m a
compressBgzfLv lv out =  do
    np <- liftIO $ getNumCapabilities
    compressBgzf' (CompressParams lv (np+2)) out

data CompressParams = CompressParams {
        compression_level :: Int,
        queue_depth :: Int }
    deriving Show

compressChunk :: Int -> Ptr CChar -> CUInt -> IO Bytes
compressChunk lv ptr len =
    allocaBytes (#{const sizeof(z_stream)}) $ \stream -> do
    buf <- mallocBytes 65536

    -- steal header from the EOF marker (length is wrong for now)
    S.unsafeUseAsCString bgzfEofMarker $ \eof ->
        forM_ [0,4..16] $ \o -> do x <- peekByteOff eof o
                                   pokeByteOff buf o (x::Word32)

    -- set up ZStream
    #{poke z_stream, msg}       stream nullPtr
    #{poke z_stream, zalloc}    stream nullPtr
    #{poke z_stream, zfree}     stream nullPtr
    #{poke z_stream, opaque}    stream nullPtr
    #{poke z_stream, next_in}   stream ptr
    #{poke z_stream, next_out}  stream (buf `plusPtr` 18)
    #{poke z_stream, avail_in}  stream len
    #{poke z_stream, avail_out} stream (65536-18-8 :: CUInt)

    z_check "deflateInit2" =<< c_deflateInit2 stream (fromIntegral lv) #{const Z_DEFLATED}
                                              (-15) 8 #{const Z_DEFAULT_STRATEGY}
    -- z_check "deflate" =<< c_deflate stream #{const Z_NO_FLUSH}
    z_check "deflate" =<< c_deflate stream #{const Z_FINISH}
    z_check "deflateEnd" =<< c_deflateEnd stream

    crc0 <- c_crc32 0 nullPtr 0
    crc  <- c_crc32 crc0 ptr len

    compressed_length <- (+) (18+8) `fmap` #{peek z_stream, total_out} stream
    when (compressed_length > 65536) $ error "produced too big a block"

    -- set length in header
    pokeByteOff buf 16 (fromIntegral $ (compressed_length-1) .&. 0xff :: Word8)
    pokeByteOff buf 17 (fromIntegral $ (compressed_length-1) `shiftR` 8 :: Word8)

    pokeByteOff buf (compressed_length-8) (fromIntegral crc :: Word32)
    pokeByteOff buf (compressed_length-4) (fromIntegral len :: Word32)

    S.unsafePackCStringFinalizer buf compressed_length (free buf)

