-- BCL to BAM command line driver
-- This should become part of Jivebunny; we'll see if the standalone
-- version will retain any value.
--
-- Basic idea:  Input is a directory of BCL files, a directory of
-- pos/locs/clocs files, and read length definitions.  Output is BAM.
--
-- If we start from a run folder, we can get the read length definitions
-- from the RunInfo.xml file, and the directories are implied.  We
-- restrict to a subset of lanes in this case.  The list of tiles can be
-- overridden from the command line.

import Bio.Bam.Header
import Bio.Bam.Trim
import Bio.Iteratee
import Bio.Iteratee.Builder
import Bio.Prelude
import Control.Exception                   ( IOException )
import Foreign.Marshal.Utils               ( copyBytes )
import Paths_jivebunny                     ( version )
import System.Console.GetOpt
import System.Directory
import System.FilePath
import System.IO                           ( hPutStrLn, stderr, stdout, withFile, IOMode(..) )
import Text.XML.Light               hiding ( Text )

import qualified Data.ByteString              as B
import qualified Data.ByteString.Unsafe       as B ( unsafeUseAsCString )
import qualified Data.Set                     as S
import qualified Data.Vector.Generic          as V
import qualified Data.Vector.Storable         as W
import qualified Data.Vector.Storable.Mutable as WM

import BCL
import Locs

-- conversion of BCL/LOCS to BAM.  We go tile by tile, and for each tile
-- we need a bunch of BCLs (one per cycle) and a LOCS.  (XXX Note on
-- parallelization:  we could read the many files in parallel, and we
-- probably should.  We could also read the files for one tile while
-- outputting another.)

tileToBam :: LaneDef -> Tile -> [ Endo BgzfTokens ]
tileToBam LaneDef{..} Tile{ tile_locs = Locs vlocs, tile_filter = Filt vfilt, ..}
    -- XXX  For the time being, require the set of cycles to be a
    -- contiguous range starting at 1.
    | tile_cycles == [ 1 .. maximum tile_cycles ]
        = map Endo $ zipWith one_cluster [0..] (V.toList vlocs)
    | otherwise = error $ "unhandled set of cycles: " ++ show tile_cycles
  where
    one_cluster i (px,py) =
        TkSetMark .
        TkWord32 maxBound .                                         -- rname
        TkWord32 maxBound .                                         -- pos
        TkWord32 (0x12480000 .|. fromIntegral lqname) .             -- lqname, mapq, bin
        TkWord32 (shiftL (flagsReadOne .|. get_flag i) 16) .        -- n_cigar, flag
        TkWord32 n_seq .                                            -- n_seq
        TkWord32 maxBound .                                         -- mrnm
        TkWord32 maxBound .                                         -- mpos
        TkWord32 0 .                                                -- isize
        qname .                                                     -- qname
        ( let Just (u,v) = cycles_read_one
          in TkBclSpecial (BclArgs BclNucsBin  tile_bcls tile_stride (u-1) (v-1) i) .
             TkBclSpecial (BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i) ) .
        indexRead "XIZ" "YIZ" cycles_index_one .
        indexRead "XJZ" "YJZ" cycles_index_two .
        TkEndRecord .

        case cycles_read_two of
            Nothing    -> id
            Just (u,v) ->
                TkSetMark .
                TkWord32 maxBound .                                         -- rname
                TkWord32 maxBound .                                         -- pos
                TkWord32 (0x12480000 .|. fromIntegral lqname) .             -- lqname, mapq, bin
                TkWord32 (shiftL (flagsReadTwo .|. get_flag i) 16) .        -- n_cigar, flag
                TkWord32 n_seq .                                            -- n_seq
                TkWord32 maxBound .                                         -- mrnm
                TkWord32 maxBound .                                         -- mpos
                TkWord32 0 .                                                -- isize
                qname .                                                     -- qname
                TkBclSpecial (BclArgs BclNucsBin  tile_bcls tile_stride (u-1) (v-1) i) .
                TkBclSpecial (BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i) .
                indexRead "XIZ" "YIZ" cycles_index_one .
                indexRead2 "XJZ" "YJZ" cycles_index_two revcom_index_two .
                TkEndRecord
      where
        lqname = B.length experiment + 5 + ilen lane_number + ilen tile_nbr + ilen px + ilen py
        qname  = TkString experiment . TkDecimal lane_number . TkDecimal tile_nbr .
                 TkDecimal (fromIntegral px) . TkDecimal (fromIntegral py) . TkWord8 0

        ilen x | x < 10 = 1
               | x < 100 = 2
               | x < 1000 = 3
               | x < 10000 = 4
               | x < 100000 = 5
               | x < 1000000 = 6
               | otherwise    = 7

        indexRead  _  _   Nothing    = id
        indexRead k1 k2 (Just (u,v)) =
              TkString k1 . TkBclSpecial (BclArgs BclNucsAsc  tile_bcls tile_stride (u-1) (v-1) i) . TkWord8 0 .
              TkString k2 . TkBclSpecial (BclArgs BclQualsAsc tile_bcls tile_stride (u-1) (v-1) i) . TkWord8 0

        indexRead2  _  _   Nothing    _ = id
        indexRead2 k1 k2 (Just (u,v)) False = indexRead k1 k2 (Just (u,v))
        indexRead2 k1 k2 (Just (u,v)) True  =
              TkString k1 . TkBclSpecial (BclArgs BclNucsAscRev  tile_bcls tile_stride (u-1) (v-1) i) . TkWord8 0 .
              TkString k2 . TkBclSpecial (BclArgs BclQualsAscRev tile_bcls tile_stride (u-1) (v-1) i) . TkWord8 0

    !n_seq = let Just (ra,re) = cycles_read_one in fromIntegral $ re-ra+1

    !flagsSingle  = fromIntegral $ flagUnmapped
    !flagsReadTwo = fromIntegral $ flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagSecondMate

    !flagsReadOne = fromIntegral $ case cycles_read_two of
            Just  _ -> flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagFirstMate
            Nothing -> flagsSingle

    get_flag j = fromIntegral $ case vfilt V.!? j of Just x | odd x -> flagFailsQC ; _ -> 0

    -- flagPaired       = 0x1
    -- flagUnmapped     = 0x4
    -- flagMateUnmapped = 0x8
    -- flagFirstMate    = 0x40
    -- flagSecondMate   = 0x80
    -- flagFailsQC      = 0x200


-- This is a somewhat hairier replacement for tileToBam.  It saves us
-- something like 10..33% wall clock time globally... That's probably
-- worth the trouble.
tileToBamBrute :: LaneDef -> Tile -> [ Endo BgzfTokens ]
tileToBamBrute LaneDef{..} Tile{ tile_locs = Locs vlocs, tile_filter = Filt vfilt, ..} =
    map go [0..V.length vlocs-1]
  where
    go i = Endo $ TkLowLevel minsize $ \bb ->
           withForeignPtr (buffer bb) $ \pp -> do
           let Just (u,v) = cycles_read_one
           rd1 <- do vv <- WM.unsafeNew (v-u+1)
                     _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                BclArgs BclNucsWide tile_bcls tile_stride (u-1) (v-1) i
                     W.unsafeFreeze vv
           qs1 <- do vv <- WM.unsafeNew (v-u+1)
                     _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i
                     W.unsafeFreeze vv

           pz <- case cycles_read_two of
                Nothing ->
                    (case find_trim default_fwd_adapters rd1 qs1 of
                        (mlen, qual1, qual2)
                            | u > v                       -> r1  0
                            | mlen == 0 && qual1 >= highq -> return
                            | qual1 < lowq                -> r1  0
                            | qual1 >= highq              -> r1t 0
                            | otherwise                   -> r1 eflagAlternative >=> r1t eflagAlternative
                          where
                            r1  = one_read flagsReadOne (u,v)
                            r1t = one_read flagsReadOne (u,u+mlen-1)
                    ) (pp `plusPtr` used bb)

                Just (u',v') -> do
                    rd2 <- do vv <- WM.unsafeNew (v'-u'+1)
                              _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                         BclArgs BclNucsWide tile_bcls tile_stride (u'-1) (v'-1) i
                              W.unsafeFreeze vv
                    qs2 <- do vv <- WM.unsafeNew (v'-u'+1)
                              _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                         BclArgs BclQualsBin tile_bcls tile_stride (u'-1) (v'-1) i
                              W.unsafeFreeze vv
                    (case find_merge default_fwd_adapters default_rev_adapters rd1 qs1 rd2 qs2 of
                        (mlen, qual1, qual2)
                             | u > v && u' > v'                     -> r1 0 >=> r2 0
                             | qual1 < lowq                         -> r1 0 >=> r2 0
                             | qual1 >= highq && mlen == 0          -> return
                             | qual1 >= highq                       -> rm 0
                             | mlen < v-u-21 || mlen < v'-u'-21     -> rm 0
                             | otherwise -> r1 eflagAlternative >=> r2 eflagAlternative >=> rm eflagAlternative

                          where
                            r1 = one_read flagsReadOne (u,v)
                            r2 = one_read flagsReadTwo (u',v')
                            rm = undefined -- XXX
                        ) (pp `plusPtr` used bb)

           if pz `minusPtr` pp - used bb <= minsize
                then return bb { used = pz `minusPtr` pp }
                else error "Oh noes, minsize was too optimistic!"
      where
        one_read :: Word32 -> (Int,Int) -> Int -> Ptr Word8 -> IO (Ptr Word8)
        one_read flags (u,v) ff p = do
            pokeByteOff p  4 (maxBound :: Word32)                               -- rname
            pokeByteOff p  8 (maxBound :: Word32)                               -- pos
            pokeByteOff p 12 (0x12480000 .|. fromIntegral lqname :: Word32)     -- lqname, mapq, bin
            pokeByteOff p 16 (shiftL (flags .|. get_flag i) 16 :: Word32)       -- n_cigar, flag
            pokeByteOff p 20 (fromIntegral (v-u+1) :: Word32)                   -- n_seq
            pokeByteOff p 24 (maxBound :: Word32)                               -- mrnm
            pokeByteOff p 28 (maxBound :: Word32)                               -- mpos
            pokeByteOff p 32 (0 :: Word32)                                      -- isize
            let ln = B.length experiment
            B.unsafeUseAsCString experiment $ \q ->
                copyBytes (p `plusPtr` 36) (castPtr q) ln
            let p0 = p `plusPtr` (36+ln)
            p1 <- plusPtr p0 <$> int_loop p0 lane_number
            p2 <- plusPtr p1 <$> int_loop p1 tile_nbr
            p3 <- plusPtr p2 <$> int_loop p2 (fromIntegral px)
            p4 <- plusPtr p3 <$> int_loop p3 (fromIntegral py)
            pokeByteOff p4 0 (0 :: Word8)

            let p5 = plusPtr p4 1
            p6 <- plusPtr p5 <$> loop_bcl_special p5 (BclArgs BclNucsBin  tile_bcls tile_stride (u-1) (v-1) i)
            p7 <- plusPtr p6 <$> loop_bcl_special p6 (BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i)

            p8 <- indexRead BclNucsAsc BclQualsAsc p7 "XIZ" "YIZ" cycles_index_one
            p9 <- if revcom_index_two
                  then indexRead BclNucsAscRev BclQualsAscRev p8 "XJZ" "YJZ" cycles_index_two
                  else indexRead BclNucsAsc    BclQualsAsc    p8 "XJZ" "YJZ" cycles_index_two
            p10 <- if ff /= 0 then do pokeByteOff p9 0 (c2w 'F')
                                      pokeByteOff p9 1 (c2w 'F')
                                      pokeByteOff p9 2 (c2w 'C')
                                      pokeByteOff p9 3 (fromIntegral ff :: Word8)
                                      return $ p9 `plusPtr` 4
                              else return p9
            pokeByteOff p 0 (fromIntegral $ (p10 `minusPtr` p) - 4 :: Word32)
            return p10


        (px,py) = vlocs V.! i
        lqname = B.length experiment + 5 + ilen lane_number + ilen tile_nbr + ilen px + ilen py

        ilen x | x < 10 = 1
               | x < 100 = 2
               | x < 1000 = 3
               | x < 10000 = 4
               | x < 100000 = 5
               | x < 1000000 = 6
               | otherwise    = 7

        indexRead :: BclSpecialType -> BclSpecialType -> Ptr Word8 -> Bytes -> Bytes -> Maybe (Int,Int) -> IO (Ptr Word8)
        indexRead  _   _  p'  _  _   Nothing    = return p'
        indexRead tp1 tp2 p' k1 k2 (Just (u,v)) = do
            B.unsafeUseAsCString k1 $ \q -> copyBytes p' (castPtr q) (B.length k1)
            l1 <- loop_bcl_special (p' `plusPtr` B.length k1) (BclArgs  tp1 tile_bcls tile_stride (u-1) (v-1) i)
            let pp1 = p' `plusPtr` (B.length k1 + l1 + 1)
            pokeByteOff pp1 (-1) (0 :: Word8)

            B.unsafeUseAsCString k2 $ \q -> copyBytes pp1 (castPtr q) (B.length k2)
            l2 <- loop_bcl_special (pp1 `plusPtr` B.length k2) (BclArgs tp2 tile_bcls tile_stride (u-1) (v-1) i)
            let pp2 = pp1 `plusPtr` (B.length k2 + l2 + 1)
            pokeByteOff pp2 (-1) (0 :: Word8)
            return pp2


    !minsize = 2 * length tile_cycles + 2 * B.length experiment + 256

    !flagsSingle  = fromIntegral $ flagUnmapped
    !flagsReadTwo = fromIntegral $ flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagSecondMate

    !flagsReadOne = fromIntegral $ case cycles_read_two of
            Just  _ -> flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagFirstMate
            Nothing -> flagsSingle

    get_flag j = fromIntegral $ case vfilt V.!? j of Just x | odd x -> flagFailsQC ; _ -> 0

    -- flagPaired       = 0x1
    -- flagUnmapped     = 0x4
    -- flagMateUnmapped = 0x8
    -- flagFirstMate    = 0x40
    -- flagSecondMate   = 0x80
    -- flagFailsQC      = 0x200

    lowq  =  20 -- XXX
    highq = 200 -- XXX


-- | Definition of a lane to be processed.  Includes paths, read
-- definitions.
data LaneDef = LaneDef {
    experiment :: !B.ByteString,
    lane_number :: !Int,

    -- | Root of BCL hierarchy, contains BCLs in subdirectories, filter
    -- and control files.
    path_bcl :: FilePath,

    -- | Path to location files.  Contains clocs, locs, or pos_txt.
    path_locs :: FilePath,

    -- | Cycles in the first business read.
    cycles_read_one :: Maybe (Int,Int),

    -- | Cycles in the second business read, if present.
    cycles_read_two :: Maybe (Int,Int),

    -- | Cycles in the first index read, if present.
    cycles_index_one :: Maybe (Int,Int),

    -- | Cycles in the second index read, if present.
    cycles_index_two :: Maybe (Int,Int),

    -- | Shall I revcom the second index sequence?
    revcom_index_two :: Bool,

    -- | List of basenames, one for each tile.
    tiles :: Maybe [String] }
  deriving Show

default_lanedef :: Int -> LaneDef
default_lanedef ln = LaneDef
    { experiment       = ""
    , lane_number      = ln
    , path_bcl         = error "need path to BCL files"
    , path_locs        = error "need path to LOCS files"
    , cycles_read_one  = Nothing
    , cycles_read_two  = Nothing
    , cycles_index_one = Nothing
    , cycles_index_two = Nothing
    , revcom_index_two = unknown_revcom_status
    , tiles            = Nothing }

unknown_revcom_status :: Bool
unknown_revcom_status = error "can't tell if index two is reverse-complemented"

data Cfg = Cfg
        { cfg_output :: (Handle -> IO ()) -> IO ()
        , cfg_report :: String -> IO ()
        -- | only used when no run folder is speficied
        , cfg_lanes :: [Int]
        -- | applied to the LaneDefs derived from a RunInfo.xml
        , cfg_overrides :: [LaneDef] -> [LaneDef] }

default_cfg :: Cfg
default_cfg = Cfg { cfg_output    = \k -> k stdout
                  , cfg_report    = const $ return ()
                  , cfg_lanes     = [1]
                  , cfg_overrides = id }

options :: [OptDescr (Cfg -> IO Cfg)]
options = [
    Option "o" ["output"]               (ReqArg set_output    "FILE") "Write output to FILE",
    Option "l" ["lanes"]                (ReqArg set_lanes     "LIST") "Process only lanes in LIST",
    Option "t" ["tiles"]                (ReqArg set_tiles     "LIST") "Process only tiles in LIST",
    Option "e" ["experiment-name"]      (ReqArg set_expname   "NAME") "Override experiment name to NAME",
    Option "b" ["bcl-path"]             (ReqArg set_bcl_path  "PATH") "Override path to BCL files",
    Option "p" ["pos-path","locs-path"] (ReqArg set_locs_path "PATH") "Override path to POS files",
    Option "r" ["read1"]                (ReqArg set_read1    "RANGE") "Read 1 comprises cycles in RANGE",
    Option "R" ["read2"]                (ReqArg set_read2    "RANGE") "Read 2 comprises cycles in RANGE",
    Option "i" ["index1"]               (ReqArg set_index1   "RANGE") "Index 1 comprises cycles in RANGE",
    Option "I" ["index2"]               (ReqArg set_index2   "RANGE") "Index 2 comprises cycles in RANGE",
    Option [ ] ["mpi-protocol"]         (NoArg         set_norevcom2) "Do not reverse-complement index read two",
    Option [ ] ["nextera-protocol"]     (NoArg           set_revcom2) "Reverse-complement index read two",
    Option "v" ["verbose"]              (NoArg           set_verbose) "Enable progress reporting",
    Option "V" ["version"]              (NoArg          disp_version) "Display program version and exit",
    Option "h?"["help","usage"]         (NoArg            disp_usage) "Display this usage information and exit" ]

  where
    set_lanes     a = override . filter $ \l -> lane_number l `elem` readWith pint_list a
    set_bcl_path  a = override $ map (\l -> l { path_bcl  = a }) . take 1
    set_locs_path a = override $ map (\l -> l { path_locs = a }) . take 1

    set_expname   a = override . map $ \l -> l { experiment       =              fromString a }
    set_read1     a = override . map $ \l -> l { cycles_read_one  = Just $ readWith  prange a }
    set_read2     a = override . map $ \l -> l { cycles_read_two  =        readWith pmrange a }
    set_index1    a = override . map $ \l -> l { cycles_index_one =        readWith pmrange a }
    set_index2    a = override . map $ \l -> l { cycles_index_two =        readWith pmrange a }

    set_revcom2     = override . map $ \l -> l { revcom_index_two =  True }
    set_norevcom2   = override . map $ \l -> l { revcom_index_two = False }

    set_output  a c = return $ c { cfg_output = \k -> withFile (a++"#") WriteMode k >> renameFile (a++"#") a }
    set_verbose   c = return $ c { cfg_report = hPutStrLn stderr }

    set_tiles   a c = override (map (\l -> l { tiles = Just . snub $ readWith pstring_list a })) $
                      c { cfg_lanes = case complete $ pint_list a of [ts] -> snub ts ; _ -> cfg_lanes c }

    disp_version _ = do pn <- getProgName
                        hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
                        exitSuccess

    disp_usage   _ = do p <- getProgName
                        hPutStrLn stderr $ "Usage: " ++ usageInfo p options
                        exitSuccess

    override f c = return $ c { cfg_overrides = f . cfg_overrides c }

snub :: Ord a => [a] -> [a]
snub = S.toList . S.fromList

complete :: [(a,String)] -> [a]
complete = map fst . filter (all isSpace . snd)

readWith :: ReadS a -> String -> a
readWith r s = case complete $ r s of
    [a] -> a ; _ -> error $ "couldn't parse " ++ show s

pmrange :: ReadS (Maybe (Int,Int))
pmrange "-" = [ (Nothing, "") ]
pmrange  s  = [ (Just r, s') | (r,s') <- prange s ]

prange :: ReadS (Int,Int)
prange s = [ ((a,b), s'') | (a,c:s') <- reads s,  c == ',' || c == '-'
                          , (b, s'') <- reads s', all isSpace s'', b >= a ]
        ++ [ ((a,a), s') | (a,s') <- reads s ]

pint_list :: ReadS [Int]
pint_list s = [ (as,s') | (as,s') <- pprim_list s ]
           ++ [ (as++as',s') | (as,',':s1) <- pprim_list s
                             , (as',s') <- pint_list s1 ]

pstring_list :: ReadS [String]
pstring_list = \s -> [ (as,s') | (as,s') <- patom s ] ++
                     [ (as++as',s') | (as,',':s1) <- patom s
                                    , (as',s') <- pstring_list s1 ]
  where
    patom = \s -> case complete $ pprim_list s of
        [ ] -> [ ([a],s') ] where (a,s') = break (==',') s
        pps -> [ (map show as,[]) | as <- pps ]

pprim_list :: ReadS [Int]
pprim_list s = [ ([a],s') | (a,s') <- reads s ]
            ++ [ ([a..b],s') | (a,'-':s1) <- reads s
                             , (b,s') <- reads s1
                             , b >= a ]

-- | Takes a run folder and derives the lane definitions.
lanesFromRun :: FilePath -> IO [LaneDef]
lanesFromRun rundir = fmap catMaybes . forM [1..8] $ \lane_number -> do
    -- Try and drop the date from the run directory.  Ignores the
    -- trailing slash if it was specified.
    let experiment = case break (== '_') (takeBaseName $ dropTrailingPathSeparator rundir) of
                        (l,r) | all isDigit l -> fromString $ drop 1 r
                        _                     -> fromString $ takeBaseName rundir

    let path_bcl = rundir </> "Data/Intensities/BaseCalls/L00" ++ [intToDigit lane_number]
        path_locs = rundir </> "Data/Intensities/L00" ++ [intToDigit lane_number]

    has_both <- (&&) <$> doesDirectoryExist path_bcl
                     <*> doesDirectoryExist path_locs

    if has_both then do
        ts     <- listTiles path_locs
        cycles <- listCycles path_bcl
        if null cycles || null ts
          then return Nothing
          else do ri <- expand_ri (maximum cycles) 1 <$> getRunInfo rundir
                  let (cycles_read_one, cycles_read_two)
                        = case map fst $ filter ((Just True /=) . snd) ri of
                            r1:r2:_ -> (Just r1, Just r2)
                            r1:_    -> (Just r1, Nothing)
                            _       -> error "shouldn't happen"

                  let (cycles_index_one, cycles_index_two)
                        = case map fst $ filter ((Just False /=) . snd) ri of
                            r1:r2:_ -> (Just r1, Just r2)
                            r1:_    -> (Just r1, Nothing)
                            _       -> (Nothing, Nothing)

                  revcom_index_two <- is_nextera_recipe rundir
                  return $ Just LaneDef{ tiles = Just ts, .. }

      else return Nothing

  where
    expand_ri total count [           ] = if count <= total then [((count,total),Nothing)] else []
    expand_ri total count ((l,isix):rs) = ((count,count+l-1),Just isix) : expand_ri total (count+l) rs

-- we deal with every tile we can find a locs file for
listTiles :: FilePath -> IO [String]
listTiles locs = snub . mapMaybe get_tile <$> getDirectoryContents locs
  where
    get_tile fn | takeExtension fn == ".gz"                  = get_tile (dropExtension fn)
                | takeExtension fn `elem` [".clocs",".locs"] = Just $ dropExtension fn
                | otherwise                                  = Nothing

listCycles :: FilePath -> IO [Int]
listCycles bcls = mapMaybe get_cycle <$> getDirectoryContents bcls
  where
    get_cycle ('C':fn) | (l,".1") <- break (== '.') fn, [(c,[])] <- reads l = Just c
    get_cycle        _                                                      = Nothing


-- For each tile, read locs and all the bcls.  Run tileToBam and emit.
bamFromBcl :: (String -> IO ()) -> LaneDef -> Enumerator (Endo BgzfTokens) IO b
bamFromBcl report ld@LaneDef{..} it0 =
    foldM (\it fn -> do report fn
                        tile <- readTile fn path_locs path_bcl [1..ce]
                        -- enumPure1Chunk (tileToBam ld tile) it
                        enumPure1Chunk (fold $ tileToBamBrute ld tile) it
          ) it0 (maybe [] id tiles)
  where
    !ce = maybe id (max . snd) cycles_index_two $
          maybe id (max . snd) cycles_index_one $
          maybe id (max . snd) cycles_read_two  $
          snd $ fromJust cycles_read_one


main :: IO ()
main = do
    (opts, rs, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Cfg{..} <- foldl (>>=) (return default_cfg) opts

    lanedefs <- case (rs, cfg_lanes) of
                    ([],[ln]) -> do let [ldef] = cfg_overrides [ default_lanedef ln ]
                                    cs <- case cycles_read_one ldef of
                                            Just zz -> return zz
                                            Nothing -> do m <- maximum <$> listCycles (path_bcl ldef)
                                                          return . (,) 1
                                                            . maybe id (max . fst) (cycles_read_two  ldef)
                                                            . maybe id (max . fst) (cycles_index_one ldef)
                                                            . maybe id (max . fst) (cycles_index_two ldef) $ m
                                    ts <- maybe (listTiles (path_locs ldef)) return $ tiles ldef
                                    return [ ldef { cycles_read_one = Just cs, tiles = Just ts } ]
                    ([],_   ) -> fail "need at least one run or exactly one lane number"
                    ( _,_   ) -> cfg_overrides . concat <$> mapM lanesFromRun rs

    mapM_ (cfg_report . show) lanedefs

    cfg_output $ \hdl ->
        (\o -> encodeHeader >>= flip enumPure1Chunk o) >=>
        foldr ((>=>) . bamFromBcl cfg_report) run lanedefs $
            joinI $ encodeBgzf 6 $ mapChunksM_ (B.hPut hdl)


-- Look for a useable XML file, either RunInfo.xml, or RunParameters.xml.
-- We'll match case insensitively, because sometimes case gets mangled
-- during the network copy.
getRunInfo :: FilePath -> IO [(Int,Bool)]
getRunInfo dir = do
    xmls <- filter (\f -> map toLower f == "runinfo.xml" || map toLower f == "runparameters.xml")
            <$> getDirectoryContents dir
    case xmls of
        fp:_ -> map snd . sort .
                mapMaybe toReadDef .
                concatMap (findChildren (unqual "Read")) .
                concatMap (findElements (unqual "Reads")) .
                maybeToList . parseXMLDoc <$> B.readFile (dir </> fp)
        [  ] -> return []

toReadDef :: Element -> Maybe (Int, (Int, Bool))
toReadDef elt = do
    nbr <- readMb =<< findAttr (unqual "Number") elt
    ncc <- readMb =<< findAttr (unqual "NumCycles") elt
    let isx = rbool $ findAttr (unqual "IsIndexedRead") elt
    return (nbr, (ncc, isx))
  where
    readMb s =  case reads (dropWhile isSpace s) of
                            [(a,b)] | all isSpace b -> Just a ; _ -> Nothing

    rbool (Just "Y") = True
    rbool (Just "y") = True
    rbool          _ = False

-- We have to allow for index2 to be reverse-complemented (Nextera
-- style) or not (MPI style).  There is obviously an option to override
-- it, and a way to detect it for a given run folder:  Check for some
-- Recipe/*.xml that doesn't match Recipe/*RunState*.  If this recipe
-- has
-- Recipe/Protocol/ChemistryRef[ChemistryName="Index2FirstBaseDark"], it
-- is a Nextera recipe.  If it has
-- Recipe/Protocol/ChemistryRef[ChemistryName="Index2Preparation"], it
-- is an MPI recipe.  Else we still don't know.  (We could infer it from
-- the order of the reads, but I consider this too fragile.)

is_nextera_recipe :: FilePath -> IO Bool
is_nextera_recipe rundir =
    do fs <- filter (not . ("runstate" `isInfixOf`) . map toLower) .
             filter (".xml" `isSuffixOf`) <$>
             getDirectoryContents (rundir </> "Recipe")
       case fs of
           [  ] -> return unknown_revcom_status
           fp:_ -> foldr match_chem_name unknown_revcom_status .
                   map (map toLower) .
                   mapMaybe (findAttr     (unqual "ChemistryName")) .
                   concatMap (findChildren (unqual "ChemistryRef")) .
                   concatMap (findChildren (unqual "Protocol")) .
                   concatMap (liftA2 (++) (findElements (unqual "Recipe"))
                                          (findElements (unqual "RecipeFile"))) .
                   maybeToList . parseXMLDoc <$> B.readFile (rundir </> "Recipe" </> fp)
    `catch` discard_io_exception
  where
    match_chem_name "index2firstbasedark" = const True
    match_chem_name "index2preparation"   = const False
    match_chem_name "indexpreparation-i5" = const False
    match_chem_name _                     = id

    discard_io_exception :: IOException -> IO Bool
    discard_io_exception _ = return unknown_revcom_status


encodeHeader :: IO (Endo BgzfTokens)
encodeHeader = do
    pn <- getProgName
    as <- getArgs

    let hdr = "@HD\tVN:1.0\n@PG\tID:" ++ pn ++ "\tPN:" ++ pn ++
              "\tCL:" ++ unwords as ++ "\tVN:" ++ showVersion version ++ "\n"

    return . Endo $ TkString "BAM\1" . TkSetMark               . TkString (fromString hdr)
                  . TkEndRecord      . TkWord32 0 {- n_refs -}

