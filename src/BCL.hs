{-# LANGUAGE CPP #-}
-- | Handling of Illumina BCL files.
-- We will support plain BCL, gzipped BCL and bgzf'ed BCL.  Plain BCL
-- starts with a cluster count (4 bytes, little-endian).  Base calls
-- follow with one byte per base:  bits [0..1] encode the base in the
-- order ACGT, bits 2..7 contain the quality score.
--
-- We will have to read from many files, so reading reasonably sized
-- blocks is imperative.  The typical bcl file on a MiSeq is 0.5MB, on a
-- HiSeq it's about 3MB.  We simply read them completely---this requires
-- 0.5-1GB of memory on a typical run, which shouldn't be a problem.
-- It's more if decompression is necessary, but still reasonable.

-- The BCLs come with a companion 'filter' file.  These start with three
-- header words:  zero, format version number, number of clusters.  The
-- remainder is one byte(!) per cluster, bit 0 is the filter flag.  We
-- expect one folder that contains the filter files and per-cycle
-- subfolders.

module BCL where

import Locs
import BasePrelude
import Zlib                             ( decompressGzip )
import Control.Concurrent.Async         ( async, wait )
import Data.Vector.Fusion.Util          ( Id )
import Data.Vector.Generic              ( unstream )
import Foreign.Ptr                      ( plusPtr, castPtr )
import Foreign.Marshal.Utils            ( copyBytes, fillBytes )
import System.Directory
import System.FilePath

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Lazy               as L
import qualified Data.ByteString.Lazy.Internal      as L ( ByteString(..) )
import qualified Data.ByteString.Unsafe             as B ( unsafeIndex, unsafeUseAsCString )
import qualified Data.Vector.Fusion.Stream.Monadic  as SS
import qualified Data.Vector.Fusion.Bundle.Monadic  as S
import qualified Data.Vector.Fusion.Bundle.Size     as S
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Storable.Mutable       as VSM
import qualified Data.Vector.Unboxed                as U

newtype BCL  = BCL  (U.Vector Word8)
newtype Filt = Filt (U.Vector Word8)

-- | We will process a tile at a time.  For reasons of cache
-- efficiency, we'll concatenate the BCLs and make sure the size of each
-- cycle is not a multiple of the cache line size.  Since we store
-- multiple cycles, we should probably track them.
--
data Tile = Tile
    { tile_nbr    :: {-# UNPACK #-} !Int
    , tile_cycles ::                [Int]
    , tile_stride :: {-# UNPACK #-} !Int
    , tile_locs   :: {-# UNPACK #-} !Locs
    , tile_filter :: {-# UNPACK #-} !Filt
    , tile_bcls   :: {-# UNPACK #-} !(VS.Vector Word8) }



-- | Reads a BCL file, which can be plain, or gzip'ed, or bgzf'ed.
-- We ignore the record count in the first quadword and always fill the
-- provided vector exactly.  Too much input data is silently ignored,
-- too little is padded liberally with zeroes.

readBCL :: VSM.IOVector Word8 -> FilePath -> IO ()
readBCL vec fp = store_to_vec vec . L.drop 4 .
                    decompressGzip . L.fromChunks . (:[]) =<< B.readFile fp

readFilt :: FilePath -> IO Filt
readFilt fp = Filt <$> readVec fp 12

readVec :: FilePath -> Int64 -> IO (U.Vector Word8)
readVec fp n = evaluate . vec_from_string . L.drop n .
                    decompressGzip . L.fromChunks . (:[]) =<< B.readFile fp

-- | Turns a lazy bytestring into a vector of words.  A straight
-- forward @fromList . toList@ would have done it, but this version
-- hopefully fuses.
vec_from_string :: L.ByteString -> U.Vector Word8
vec_from_string = unstream . S.concatMap stream_bs . stream_lbs
  where
    stream_bs :: B.ByteString -> S.Bundle Id v Word8
    stream_bs bs = S.fromStream (SS.Stream step 0) (S.Exact $ B.length bs)
      where
        step i | i == B.length bs = return $ SS.Done
               | otherwise        = return $ SS.Yield (B.unsafeIndex bs i) (i+1)

    stream_lbs :: L.ByteString -> S.Bundle Id v B.ByteString
    stream_lbs lbs = S.fromStream (SS.Stream step lbs) S.Unknown
      where
        step  L.Empty       = return $ SS.Done
        step (L.Chunk c cs) = return $ SS.Yield c cs

store_to_vec :: VSM.IOVector Word8 -> L.ByteString -> IO ()
store_to_vec v ls = VSM.unsafeWith v                                $ \p ->
                    cp (castPtr p) (VSM.length v) (L.toChunks ls)
  where
    cp pd l (s:ss)
        | l >= B.length s = do B.unsafeUseAsCString s $ \ps ->
                                    copyBytes pd ps (B.length s)
                               cp (plusPtr pd (B.length s)) (l - B.length s) ss

        | otherwise       =    B.unsafeUseAsCString s $ \ps ->
                                    copyBytes pd ps l

    cp pd l [] = fillBytes pd 0 l



-- | Read a subset of the cycles of a tile.  We read the locations
-- first, because that tells us the number of clusters.  We then
-- allocate a giant vector for the contents of the BCL files, after
-- rounding the cluster number up to an odd number.
-- Arguments are the name of the locs or clocs file without extension,
-- the path to the bcl and filter file, list of cycles, ...
--
-- XXX  We use a storable vector (should be fine, it's hardly more than
-- a pointer) so we can copy directly into it.  The gunzip code could
-- make use of that, but right now doesn't.

readTile :: FilePath -> FilePath -> FilePath -> [Int] -> IO Tile
readTile fn path_locs path_bcl tile_cycles =
  do
    tile_locs   <- get_locs
    tile_filter <- get_filt

    let tile_stride = num_locs tile_locs
    vec <- VSM.new (length tile_cycles * tile_stride)

    mapM_ wait <=< sequence $ zipWith
        (\off ncycle -> let fn_bcl = path_bcl ++ "/C" ++ show ncycle ++ ".1/" ++ fn ++ ".bcl"
                            vec'   = VSM.slice off (num_locs tile_locs) vec
                        in async $ try_read_or (fn_bcl <.> "gz") (readBCL vec') $
                                   try_read_or fn_bcl            (readBCL vec') $
                                   VSM.set vec' 0

        ) [0, tile_stride ..] tile_cycles

    tile_bcls <- VS.unsafeFreeze vec
    return $ Tile {..}

  where
    tile_nbr = case reads . reverse . takeWhile (/= '_')
                  . reverse $ takeDirectory fn_locs of
                        [(n,"")] -> n ; _ -> 0

    fn_locs  = path_locs </> fn
    get_locs = try_read_or (fn_locs <.> "clocs.gz") readClocs $
               try_read_or (fn_locs <.> "clocs")    readClocs $
               try_read_or (fn_locs <.> "locs.gz")  readLocs  $
               try_read_or (fn_locs <.> "locs")     readLocs  $
               return                              (Locs U.empty)

    fn_filt  = path_bcl </> path_locs <.> "filter"
    get_filt = try_read_or (fn_filt <.> "gz") readFilt $
               try_read_or  fn_filt           readFilt $
               return                        (Filt U.empty)

    try_read_or f r k = do e <- doesFileExist f
                           if e then r f else k
