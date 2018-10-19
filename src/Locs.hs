-- | Parsing Illumina location files.  It appears we have to support
-- pos.txt, locs and clocs files; don't know if they come in compressed
-- versions, too.  Either way, we have one file per tile, and it lists X
-- and Y coordinates for each cluster.  We turn them into 16-bit integers.
--
-- Unfortunately, these files don't have any easy way to recognize the
-- format.  So we expose three different readers and a higher level
-- needs to decide which one to call.

module Locs where

import Bio.Prelude
import Bio.Util.Zlib                        ( decompressGzip )
import Data.Vector.Fusion.Bundle.Monadic    ( fromStream )
import Data.Vector.Fusion.Bundle.Size       ( Size(..) )
import Data.Vector.Fusion.Stream.Monadic    ( Stream(..), Step(..) )
import Data.Vector.Generic                  ( unstream )
import Foreign.C.Types                      ( CChar, CDouble(..) )
import Foreign.Marshal.Alloc                ( alloca )
import Foreign.Marshal.Utils                ( copyBytes )

import qualified Data.ByteString                 as B
import qualified Data.ByteString.Unsafe          as B
import qualified Data.ByteString.Lazy            as L
import qualified Data.ByteString.Lazy.Char8      as C
import qualified Data.Vector.Storable            as S
import qualified Data.Vector.Storable.Mutable    as SM
import qualified Data.Vector.Unboxed             as U

newtype Locs = Locs (U.Vector (Word32,Word32))

num_locs :: Locs -> Int
num_locs (Locs v) = U.length v

-- Pos files are text, the first two words on each line are x and y
-- coordinate (signed decimal floating point numbers).  They are rounded
-- to integers and clamped to a minimum of zero.
readPosTxt :: FilePath -> IO Locs
readPosTxt fp = Locs . U.fromList . map parse . C.lines . decompressGzip <$> L.readFile fp
  where
    round' :: CDouble -> Word32
    round' = round . max 0

    parse l = let s = B.concat $ C.toChunks l ++ ["\0"]
              in unsafePerformIO $ B.unsafeUseAsCString s $ \p0  ->
                 alloca                                   $ \ppe ->
                 strtod p0 ppe                          >>= \x   ->
                 peek ppe                               >>= \p1  ->
                 strtod p1 ppe                          >>= \y   ->
                 return (round' x, round' y)

foreign import ccall unsafe "static stdlib.h strtod" strtod
    :: Ptr CChar -> Ptr (Ptr CChar) -> IO CDouble


-- | Locs files have three header words (4 bytes, little endian), the
-- third is the number of clusters.  They are followed by two floats
-- (IEEE single precision, little endian) for (x,y) of each cluster.  We
-- round them off and clamp at zero.

readLocs :: FilePath -> IO Locs
readLocs fp = fmap conv . vec_from_lbs . L.drop 12 . decompressGzip =<< L.readFile fp
  where
    conv vec = Locs $ U.unfoldrN (S.length vec) step (vec :: S.Vector Float)
    step vec | S.length vec < 2 = Nothing
             | otherwise        = Just ((round' $ vec S.! 0, round' $ vec S.! 1), S.drop 2 vec)
    round' = round . max 0 . (+) 1000 . (*) 10


-- clocs files store position data for successive clusters, compressed in bins as follows:
--     Byte 0   : unused
--     Byte 1-4 : unsigned int numBins
--     The rest of the file consists of bins/blocks, where a bin consists of an integer
--     indicating number of blocks, followed by that number of blocks and a block consists
--     of an x-y coordinate pair.  In otherwords:
--
--     for each bin
--         byte 1: Unsigned int numBlocks
--         for each block:
--             byte 1 : byte xRelativeCoordinate
--             byte 2 : byte yRelativeCoordinate
--
--     Actual x and y values are computed using the following algorithm
--
--     xOffset = yOffset = 0
--     imageWidth = 2048
--     blockSize = 25
--     maxXbins = ceil(imageWidth/blockSize)
--     for each bin:
--         for each block:
--             x = xRelativeCoordinate/10 + xoffset
--             y = yRelativeCoordinate/10 + yoffset
--
--         if (binIndex > 0 && ((binIndex + 1) % maxXbins == 0)):
--            xOffset = 0
--            yOffset += blockSize
--         else:
--            xOffset += blockSize-
--
-- (what an ugly read)

data CLocs_Args = CLocs_Args !L.ByteString !Word32 !Word32 !Word8 !Int

readClocs :: FilePath -> IO Locs
readClocs fp = Locs . unstream . (\inp -> fromStream (Stream (return . decode) (CLocs_Args inp 0 0 0 (-1))) Unknown)
                    . L.drop 5 . decompressGzip <$> L.readFile fp
  where
    imageWidth = 20480
    blockSize  = 250
    maxXbins   = fromIntegral $ (imageWidth + blockSize -1) `div` blockSize

    decode (CLocs_Args inp xo yo 0 ibin)
        = case L.uncons inp of
            Just (numb, inp')

                | succ ibin == 0
                    -> Skip $! CLocs_Args inp' 1000 1000 numb (succ ibin)

                | succ ibin `mod` maxXbins == 0
                    -> Skip $! CLocs_Args inp' 1000 (yo + blockSize) numb (succ ibin)

                | otherwise
                    -> Skip $! CLocs_Args inp' (xo + blockSize) yo numb (succ ibin)

            Nothing -> Done

    decode (CLocs_Args inp xo yo numb ibin)
        = case L.uncons inp of
            Just (xrel, inp1) -> case L.uncons inp1 of
                Just (yrel, inp2) ->
                    Yield ( xo + fromIntegral xrel, yo + fromIntegral yrel )
                          $! CLocs_Args inp2 xo yo (pred numb) ibin
                Nothing -> Done
            Nothing -> Done


vec_from_lbs :: Storable a => L.ByteString -> IO (S.Vector a)
vec_from_lbs lbs = do
    v <- SM.new $ fromIntegral $ L.length lbs
    _ <- SM.unsafeWith (v :: SM.IOVector CChar) $ \pv ->
            foldM step pv (L.toChunks lbs)
    S.unsafeFreeze (SM.unsafeCast v)
  where
    step pv str =
        B.unsafeUseAsCStringLen str $ \(ps,len) -> do
            copyBytes pv ps len
            return $! plusPtr pv len
