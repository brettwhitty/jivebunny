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
import Control.Monad.Trans.State.Strict
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

type MergeConf = (Int, Maybe Int)

-- | This is a somewhat hairier replacement for the original algorithm
-- that produced proper 'BgzfTokens'.  It saves us something like
-- 10..33% wall clock time globally... That's probably worth the
-- trouble.  (If merging is enabled, it doesn't make a measurable
-- difference, because merging is so much more expensive.)
tileToBam :: MergeConf -> LaneDef -> Tile -> [ Endo BgzfTokens ]
tileToBam merge_conf LaneDef{..} Tile{ tile_locs = Locs vlocs, tile_filter = Filt vfilt, ..}
    = map go [0..V.length vlocs-1]
  where
    go i = Endo $ TkLowLevel minsize $ \bb ->
           withForeignPtr (buffer bb) $ \pp -> do
           pz <- execStateT (go' i merge_conf) (pp `plusPtr` used bb)
           if pz `minusPtr` pp - used bb <= minsize
                then return bb { used = pz `minusPtr` pp }
                else error "Oh noes, minsize was too optimistic!"

    go' i (_, Nothing) = do
           let Just (u,v) = cycles_read_one
           one_read i Nothing flagsReadOne (u,v) 0
           forM_ cycles_read_two $ \cs -> one_read i Nothing flagsReadTwo cs 0

    go' i (low_qual, Just high_qual) = do
           let Just (u,v) = cycles_read_one
           rd1 <- lift $ do vv <- WM.unsafeNew (v-u+1)
                            _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                       BclArgs BclNucsWide tile_bcls tile_stride (u-1) (v-1) i
                            W.unsafeFreeze vv
           qs1 <- lift $ do vv <- WM.unsafeNew (v-u+1)
                            _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                       BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i
                            W.unsafeFreeze vv

           case cycles_read_two of
                Nothing ->
                    case find_trim default_fwd_adapters rd1 qs1 of
                        (mlen, qual1, qual2)
                            | u > v                            -> r1  0
                            | mlen == 0 && qual1 >= high_qual  -> return ()
                            | qual1 <  low_qual                -> r1  0
                            | qual1 >= high_qual               -> r1t eflagTrimmed
                            | otherwise -> r1 eflagAlternative >> r1t (eflagAlternative .|. eflagTrimmed)
                          where
                            r1  = one_read i (Just (qual1, qual2)) flagsReadOne (u,v)
                            r1t = one_read i (Just (qual1, qual2)) flagsReadOne (u,u+mlen-1)

                Just (u',v') -> do
                    rd2 <- lift $ do vv <- WM.unsafeNew (v'-u'+1)
                                     _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                                BclArgs BclNucsWide tile_bcls tile_stride (u'-1) (v'-1) i
                                     W.unsafeFreeze vv
                    qs2 <- lift $ do vv <- WM.unsafeNew (v'-u'+1)
                                     _  <- WM.unsafeWith vv $ \p -> loop_bcl_special (castPtr p) $
                                                BclArgs BclQualsBin tile_bcls tile_stride (u'-1) (v'-1) i
                                     W.unsafeFreeze vv
                    case find_merge default_fwd_adapters default_rev_adapters rd1 qs1 rd2 qs2 of
                        (mlen, qual1, qual2)
                             | u > v && u' > v'                     -> r1 0 >> r2 0
                             | qual1 <  low_qual                    -> r1 0 >> r2 0
                             | qual1 >= high_qual && mlen == 0      -> return ()
                             | qual1 >= high_qual                   -> rm 0
                             | mlen < v-u-21 || mlen < v'-u'-21     -> rm 0
                             | otherwise -> r1 eflagAlternative >> r2 eflagAlternative >> rm eflagAlternative

                          where
                            r1 = one_read i (Just (qual1, qual2)) flagsReadOne (u,v)
                            r2 = one_read i (Just (qual1, qual2)) flagsReadTwo (u',v')

                            rm ff = one_read_with merge_core i (Just (qual1, qual2)) flagsSingle (u,u+mlen-1) $
                                        ff .|. eflagMerged .|. (if mlen < v-u+1 then eflagTrimmed else 0)

                            merge_core _ _ _ = do
                                let qmax = fromIntegral $ min 63 qual1
                                V.imapM_  n4        $ merged_seq       mlen rd1 qs1 rd2 qs2
                                V.mapM_  (w8 . unQ) $ merged_qual qmax mlen rd1 qs1 rd2 qs2

    s8 :: B.ByteString -> StateT (Ptr Word8) IO ()
    s8 s = StateT $ \p -> B.unsafeUseAsCString s $ \q ->
                ((),p `plusPtr` B.length s) <$ copyBytes p (castPtr q) (B.length s)

    w8 :: Word8 -> StateT (Ptr Word8) IO ()
    w8 w = StateT $ \p -> ((),p `plusPtr` 1) <$ pokeByteOff p 0 w

    n4 :: Int -> Nucleotides -> StateT (Ptr Word8) IO ()
    n4 i (Ns n) | even i    = StateT $ \p -> ((),p) <$ pokeByteOff p 0 (shiftL n 4)
                | otherwise = StateT $ \p -> do a <- peekByteOff p 0
                                                pokeByteOff p 0 (a .|. n)
                                                return ((),p `plusPtr` 1)

    w32 :: Word32 -> StateT (Ptr Word8) IO ()
    w32 w = StateT $ \p -> ((),p `plusPtr` 4) <$ pokeByteOff p 0 w

    wrap :: (Ptr a -> t -> IO Int) -> t -> StateT (Ptr a) IO ()
    wrap f a = StateT $ \p -> f p a >>= \l -> return ((), p `plusPtr` l)

    one_read :: Int -> Maybe (Int,Int) -> Word32 -> (Int,Int) -> Int -> StateT (Ptr Word8) IO ()
    one_read = one_read_with $ \i u v -> do
        wrap loop_bcl_special (BclArgs BclNucsBin  tile_bcls tile_stride (u-1) (v-1) i)
        wrap loop_bcl_special (BclArgs BclQualsBin tile_bcls tile_stride (u-1) (v-1) i)

    one_read_with core i mmqs flags (u,v) ff = do
        let (px,py) = vlocs V.! i
        let lqname = B.length experiment + 5 + ilen lane_number + ilen tile_nbr + ilen px + ilen py

        p0 <- get
        modify (`plusPtr` 4)
        w32 maxBound                                                        -- rname
        w32 maxBound                                                        -- pos
        w32 (0x12480000 .|. fromIntegral lqname)                            -- lqname, mapq, bin
        w32 (shiftL (flags .|. get_flag i) 16)                              -- n_cigar, flag
        w32 (fromIntegral (v-u+1))                                          -- n_seq
        w32 maxBound                                                        -- mrnm
        w32 maxBound                                                        -- mpos
        w32 0                                                               -- isize

        s8 experiment                                                       -- rname, many pieces
        wrap int_loop lane_number
        wrap int_loop tile_nbr
        wrap int_loop (fromIntegral px)
        wrap int_loop (fromIntegral py)
        w8 0

        _ <- core i u v

        indexRead i BclNucsAsc BclQualsAsc "XIZ" "YIZ" cycles_index_one
        if revcom_index_two
              then indexRead i BclNucsAscRev BclQualsAscRev "XJZ" "YJZ" cycles_index_two
              else indexRead i BclNucsAsc    BclQualsAsc    "XJZ" "YJZ" cycles_index_two

        when (ff /= 0) $ s8 "FFC" >> w8 (fromIntegral ff)

        forM_ mmqs $ \(ym,yn) -> do s8 "YMC" >> w8 (fromIntegral (min 255 ym))
                                    s8 "YNC" >> w8 (fromIntegral (min 255 yn))
        StateT $ \pe ->
            ((),pe) <$ pokeByteOff p0 0 (fromIntegral $ (pe `minusPtr` p0) - 4 :: Word32)


    ilen x | x < 10 = 1
           | x < 100 = 2
           | x < 1000 = 3
           | x < 10000 = 4
           | x < 100000 = 5
           | x < 1000000 = 6
           | otherwise    = 7

    indexRead :: Int -> BclSpecialType -> BclSpecialType -> Bytes -> Bytes -> Maybe (Int,Int) -> StateT (Ptr Word8) IO ()
    indexRead _  _   _   _  _   Nothing    = StateT $ \p' -> return ((),p')
    indexRead i tp1 tp2 k1 k2 (Just (u,v)) = StateT $ \p' -> do
        B.unsafeUseAsCString k1 $ \q -> copyBytes p' (castPtr q) (B.length k1)
        l1 <- loop_bcl_special (p' `plusPtr` B.length k1) (BclArgs  tp1 tile_bcls tile_stride (u-1) (v-1) i)
        let pp1 = p' `plusPtr` (B.length k1 + l1 + 1)
        pokeByteOff pp1 (-1) (0 :: Word8)

        B.unsafeUseAsCString k2 $ \q -> copyBytes pp1 (castPtr q) (B.length k2)
        l2 <- loop_bcl_special (pp1 `plusPtr` B.length k2) (BclArgs tp2 tile_bcls tile_stride (u-1) (v-1) i)
        let pp2 = pp1 `plusPtr` (B.length k2 + l2 + 1)
        pokeByteOff pp2 (-1) (0 :: Word8)
        return ((),pp2)


    minsize = 4 * length tile_cycles + 3 * B.length experiment + 512

    flagsSingle  = fromIntegral $ flagUnmapped
    flagsReadTwo = fromIntegral $ flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagSecondMate

    flagsReadOne = case cycles_read_two of
            Just  _ -> fromIntegral $ flagUnmapped .|. flagMateUnmapped .|. flagPaired .|. flagFirstMate
            Nothing -> flagsSingle

    get_flag j = fromIntegral $ case vfilt V.!? j of Just x | odd x -> flagFailsQC ; _ -> 0


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
        { cfg_output    :: (Handle -> IO ()) -> IO ()
        , cfg_report    :: String -> IO ()
        -- | only used when no run folder is speficied
        , cfg_lanes     :: [Int]
        -- | applied to the LaneDefs derived from a RunInfo.xml
        , cfg_overrides :: [LaneDef] -> [LaneDef]
        , cfg_merge     :: MergeConf }

default_cfg :: Cfg
default_cfg = Cfg { cfg_output    = \k -> k stdout
                  , cfg_report    = const $ return ()
                  , cfg_lanes     = [1]
                  , cfg_overrides = id
                  , cfg_merge     = (20, Nothing) }

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
    Option "m" ["merge-overlap"]        (OptArg set_merge     "QUAL") "Attempt to merge or trim reads",
    Option "q" ["merge-qual"]           (ReqArg set_qual      "QUAL") "Minimum quality for merge is QUAL",
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

    set_merge Nothing  c =                    return $ c { cfg_merge = ( fst (cfg_merge c), Just 200 ) }
    set_merge (Just a) c = readIO a >>= \m -> return $ c { cfg_merge = ( fst (cfg_merge c), Just m ) }
    set_qual        a  c = readIO a >>= \q -> return $ c { cfg_merge = ( q, snd (cfg_merge c) ) }

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
bamFromBcl :: (String -> IO ()) -> MergeConf -> LaneDef -> Enumerator (Endo BgzfTokens) IO b
bamFromBcl report merge_conf ld@LaneDef{..} it0 =
    foldM (\it fn -> do report fn
                        tile <- readTile fn path_locs path_bcl [1..ce]
                        enumPure1Chunk (fold $ tileToBam merge_conf ld tile) it
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
        foldr ((>=>) . bamFromBcl cfg_report cfg_merge) run lanedefs $
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

