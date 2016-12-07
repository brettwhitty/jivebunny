-- Two-stage demultiplexing.
--
-- We assume we know the list of i7 and i5 index oligos.  We seek to
-- decompose a set of reads into a mix of pairs of these by the Maximum
-- Likelihood method.  Once that's done, an empirical Bayesian Maximum
-- Posterior call is done.  All kinds of errors can be rolled into one
-- quality score.
--
--  - Input layer to gather index sequences.  (Done.)
--  - Input layer to gather read group defs.  (Done.)
--  - First pass to gather data.  Any index read shall be represented
--    in a single Word64.  (Done.  Reading BAM is slow.  BCL would be
--    much more suitable here.)
--  - Multiple passes of the EM algorithm.  (Done.)
--  - Start with a naive mix, to avoid arguments.  (Done.)
--  - Final calling pass from BAM to BAM.  (Done.  BCL to BAM would be
--    even nicer.)
--  - Auxillary statistics:  composition of the mix (Done.), false
--    assignment rates per read group (Done.), maximum achievable false
--    assignment rates (Done.)

import Bio.Bam
import Bio.Prelude
import Bio.Util.Numeric                 ( showNum )
import Data.Aeson
import Foreign.C.Types
import Foreign.Marshal.Alloc
import Foreign.Ptr
import Foreign.Storable
import Paths_biohazard                  ( version, getDataFileName )
import System.Console.GetOpt
import System.Random                    ( randomRIO )

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as BS
import qualified Data.HashMap.Strict                as HM
import qualified Data.Text                          as T
import qualified Data.Text.Encoding                 as T
import qualified Data.Text.IO                       as T
import qualified Data.Text.Lazy                     as L hiding ( singleton )
import qualified Data.Text.Lazy.IO                  as L
import qualified Data.Text.Lazy.Builder             as L
import qualified Data.Text.Lazy.Builder.Int         as L
import qualified Data.Text.Lazy.Builder.RealFloat   as L
import qualified Data.Vector                        as V
import qualified Data.Vector.Algorithms.Intro       as V
import qualified Data.Vector.Unboxed                as U
import qualified Data.Vector.Storable               as VS
import qualified Data.Vector.Storable.Mutable       as VSM
import qualified Data.Vector.Generic                as VG
import qualified Data.Vector.Generic.Mutable        as VGM

import Index

fromS :: B.ByteString -> Index
fromS sq = fromSQ sq (B.replicate (B.length sq) 64)

fromSQ :: B.ByteString -> B.ByteString -> Index
fromSQ sq qs = Index . foldl' (\a b -> a `shiftL` 8 .|. fromIntegral b) 0 $
               take 8 $ (++ repeat 0) $
               B.zipWith (\b q -> shiftL (b .&. 0xE) 4 .|. (min 31 $ max 33 q - 33)) sq qs

fromTags :: BamKey -> BamKey -> BamRaw -> Index
fromTags itag qtag br = fromSQ sq  (if B.null qs then "@@@@@@@@" else qs)
  where
    sq = extAsString itag $ unpackBam br
    qs = extAsString qtag $ unpackBam br

gather :: MonadIO m => Int -> (String -> IO ()) -> (String -> IO ()) -> BamMeta -> Iteratee [BamRaw] m (U.Vector (Index, Index))
gather num say mumble hdr = case hdr_sorting $ meta_hdr hdr of
    Unsorted    -> greedy
    Grouped     -> greedy
    Queryname   -> greedy
    Unknown     -> safe
    Coordinate  -> fair
    GroupSorted -> fair
  where
    greedy = do liftIO . say $ "File is unsorted, sampling up to "
                            ++ showNum num ++ " records from the beginning.\n"
                go stream2vectorN

    fair   = do liftIO . say $ "File is sorted, need to sample up to "
                            ++ showNum num ++ " from whole file.\n"
                go subsam2vector

    safe   = do liftIO . say $ "File might be sorted, need to sample up to "
                            ++ showNum num ++ " from whole file.\n"
                go subsam2vector

    go k = filterStream ((\b -> not (isPaired b) || isFirstMate b) . unpackBam) =$
           progressNum "reading " 0x100000 mumble =$
           mapStream (fromTags "XI" "YI" &&& fromTags "XJ" "YJ") =$ k num


subsam2vector :: (MonadIO m, ListLike s a, Nullable s, VG.Vector v a) => Int -> Iteratee s m (v a)
subsam2vector sz = liftIO (VGM.new sz) >>= go 0
  where
    go !i !mv = tryHead >>= \x -> case x of
                  Nothing -> liftIO $ if i < sz then VG.unsafeFreeze $ VGM.take i mv
                                                else VG.unsafeFreeze mv
                  Just  a -> do liftIO $ if i < sz
                                    then VGM.write mv i a
                                    else do p <- randomRIO (0,i)
                                            when (p < sz) $ VGM.write mv p a
                                go (i+1) mv

data IndexTab = IndexTab { unique_indices :: U.Vector Index
                         , canonical_names :: V.Vector T.Text
                         , alias_names :: HM.HashMap T.Text Int }

single_placeholder :: IndexTab
single_placeholder = IndexTab (U.singleton (fromS "")) (V.singleton "is4") $
                        HM.fromList [ (k,0) | [_,_,k] <- map T.words $ T.lines default_rgs ]

data Both = Both { p7is :: IndexTab, p5is :: IndexTab }

instance FromJSON Both where
    parseJSON = withObject "toplevel object expected" $ \v ->
                          both <$>  ((v .: "p7index") >>= parse_assocs)
                               <*> (((v .: "p5index") >>= parse_assocs) <|> pure [])
      where
        parse_assocs = withObject "association list expected" $ \o ->
                            sequence [ (,) k <$> withText "sequence expected" (return . T.encodeUtf8) v | (k,v) <- HM.toList o ]

        both as7 as5 = Both (canonical as7) (canonical as5)
          where
            canonical pps =
                let hm = HM.toList $ HM.fromListWith (++) [ (fromS v,[k]) | (k,v) <- pps ]
                in IndexTab (U.fromList $ map fst hm)
                            (V.fromList $ map (head . snd) hm)
                            (HM.fromList $ [ (k,i) | (i, ks) <- zip [0..] (map snd hm), k <- ks ])

data RG = RG { rgid :: B.ByteString
             , rgi7 :: Int
             , rgi5 :: Int
             , tags :: BamOtherShit }

-- | Parses read group defintions from a file.  The file can have
-- optional header lines, the remainder must be a tab-separated table,
-- first column is the read group name, second is the P7 index name,
-- third is the P5 index name (*must* be present), all others are tagged
-- fields just like BAM expects them in the header.
--
-- For integration with a LIMS, something structured like JSON would
-- probably work better, however, absent such a LIMS, tables are easier
-- to come by.

readRGdefns :: HM.HashMap T.Text Int -> HM.HashMap T.Text Int -> T.Text -> [ RG ]
readRGdefns p7is p5is = map repack . filter (not . null) . map (T.split (=='\t'))
                      . dropWhile ("#" `T.isPrefixOf`) . T.lines
  where
    repack (rg:_) | T.any (\c -> c == '/' || c == ',') rg = error $ "RG name must not contain ',' or '/': " ++ show rg
    repack (rg:p7:p5:tags) = case HM.lookup p7 p7is of
        Nothing -> error $ "unknown P7 index " ++ show p7
        Just i7 -> case HM.lookup p5 p5is of
            Nothing -> error $ "unknown P5 index " ++ show p5
            Just i5 -> RG (T.encodeUtf8 rg) i7 i5 (map repack1 tags)
    repack ws = error $ "short RG line " ++ show (T.intercalate "\t" ws)
    repack1 w | T.length w > 3 && T.index w 2 == ':'
                    = (fromString [T.index w 0, T.index w 1], T.encodeUtf8 $ T.drop 3 w)
              | otherwise = error $ "illegal tag " ++ show w

default_rgs :: T.Text
default_rgs = "PhiXA\tPhiA\tPhiA\nPhiXC\tPhiC\tPhiC\nPhiXG\tPhiG\tPhiG\nPhiXT\tPhiT\tPhiT\nPhiX\tPhiX\tis4\n"

-- | Compute mismatch score: sum of the qualities in 'a' at positions
-- where the bases don't match.  Works by comparing through an xor,
-- building a mask from it, then adding quality scores sideways.
--
-- Since we keep quality scores in the lower 5 bits of each byte, adding
-- all eight is guaranteed to fit into the highest 8 bits.
match :: Index -> Index -> Word64
match (Index a) (Index b) = score
  where x = a `xor` b
        y = (shiftR x 5 .|. shiftR x 6 .|. shiftR x 7) .&. 0x0101010101010101
        bitmask = (0x2020202020202020 - y) .&. 0x1F1F1F1F1F1F1F1F
        score = shiftR ((a .&. bitmask) * 0x0101010101010101) 56

-- | A mixture description is one probability for each combination of p7
-- and p5 index.  They should sum to one.
type Mix = VS.Vector Double
type MMix = VSM.IOVector Double

padding :: Int
padding = 31

stride' :: Int -> Int
stride' n5 = (n5 + padding) .&. complement padding

-- | Computing the naively assumed mix when nothing is known:  uniform
-- distribution.
naiveMix :: (Int,Int) -> Int -> Mix
naiveMix (n7,n5) total = VS.replicate vecsize (fromIntegral total / fromIntegral bins)
  where
    !vecsize = n7 * stride' n5
    !bins    = n7 * n5

-- | Matches an index against both p7 and p5 lists, computes posterior
-- likelihoods from the provided prior and accumulates them onto the
-- posterior.
unmix1 :: U.Vector Index -> U.Vector Index -> Mix -> MMix -> (Index, Index) -> IO ()
unmix1 p7 p5 prior acc (x,y) =
    let !m7 = VS.fromListN (U.length p7) . map (phredPow . match x) $ U.toList p7
        !l5 = stride' (U.length p5)
        !m5 = VS.fromListN l5 $ map (phredPow . match y) (U.toList p5) ++ repeat 0

    -- *sigh*, Vector doesn't fuse well.  Gotta hand it over to gcc.  :-(
    in VSM.unsafeWith acc                                           $ \pw ->
       VS.unsafeWith prior                                          $ \pv ->
       VS.unsafeWith m7                                             $ \q7 ->
       VS.unsafeWith m5                                             $ \q5 ->
       c_unmix_total pv q7 (fromIntegral $ VS.length m7)
                        q5 (fromIntegral $ l5 `div` succ padding)
                        nullPtr nullPtr                           >>= \total ->
       c_unmix_qual pw pv q7 (fromIntegral $ VS.length m7)
                          q5 (fromIntegral $ l5 `div` succ padding)
                          total 0 0                               >>= \_qual ->
       return ()    -- the quality is meaningless here

foreign import ccall unsafe "c_unmix_total"
    c_unmix_total :: Ptr Double                     -- prior mix
                  -> Ptr Double -> CUInt            -- P7 scores, length
                  -> Ptr Double -> CUInt            -- P5 scores, length/32
                  -> Ptr CUInt -> Ptr CUInt         -- out: ML P7 index, P5 index
                  -> IO Double                      -- total likelihood

foreign import ccall unsafe "c_unmix_qual"
    c_unmix_qual :: Ptr Double                      -- posterior mix, mutable accumulator
                 -> Ptr Double                      -- prior mix
                 -> Ptr Double -> CUInt             -- P7 scores, length
                 -> Ptr Double -> CUInt             -- P5 scores, length/32
                 -> Double                          -- total likelihood
                 -> CUInt -> CUInt                  -- maximizing P7 index, P5 index
                 -> IO Double                       -- posterior probability for any other assignment

-- | Matches an index against both p7 and p5 lists, computes MAP
-- assignment and quality score.
class1 :: HM.HashMap (Int,Int) (B.ByteString, VSM.IOVector Double)
       -> U.Vector Index -> U.Vector Index
       -> Mix -> (Index, Index) -> IO (Double, Int, Int)
class1 rgs p7 p5 prior (x,y) =
    let !m7 = VS.fromListN (U.length p7) . map (phredPow . match x) $ U.toList p7
        !l5 = stride' (U.length p5)
        !m5 = VS.fromListN l5 $ map (phredPow . match y) (U.toList p5) ++ repeat 0

    -- *sigh*, Vector doesn't fuse well.  Gotta hand it over to gcc.  :-(
    in alloca                                                       $ \pi7 ->
       alloca                                                       $ \pi5 ->
       VS.unsafeWith prior                                          $ \pv ->
       VS.unsafeWith m7                                             $ \q7 ->
       VS.unsafeWith m5                                             $ \q5 ->
       ( {-# SCC "c_unmix_total" #-}
         c_unmix_total pv q7 (fromIntegral $ VS.length m7)
                          q5 (fromIntegral $ l5 `div` succ padding)
                          pi7 pi5 )                               >>= \total ->
       peek pi7                                                   >>= \i7 ->
       peek pi5                                                   >>= \i5 ->
       withDirt (fromIntegral i7, fromIntegral i5)                  $ \pw ->
       ( {-# SCC "c_unmix_qual" #-}
         c_unmix_qual pw pv q7 (fromIntegral $ VS.length m7)
                            q5 (fromIntegral $ l5 `div` succ padding)
                            total i7 i5 )                         >>= \qual ->
       return ( qual, fromIntegral i7, fromIntegral i5 )
  where
    withDirt ix k = case HM.lookup ix rgs of
            Just (_,dirt) -> VSM.unsafeWith dirt k
            Nothing       -> k nullPtr


phredPow :: Word64 -> Double
phredPow x = exp $ -0.1 * log 10 * fromIntegral x

-- | One iteration of the EM algorithm.  Input is a vector of pairs of
-- indices, the p7 and p5 index collections, and a prior mixture; output
-- is the posterior mixture.
iterEM :: U.Vector (Index, Index) -> U.Vector Index -> U.Vector Index -> Mix -> IO Mix
iterEM pps p7 p5 prior = do
    acc <- VSM.replicate (VS.length prior) 0
    U.mapM_ (unmix1 p7 p5 prior acc) pps
    VS.unsafeFreeze acc

data Loudness = Quiet | Normal | Loud

unlessQuiet :: Monad m => Loudness -> m () -> m ()
unlessQuiet Quiet _ = return ()
unlessQuiet     _ k = k

-- should I have a config for merging here?  adapter lists?
-- does it ever make sense to skip the merging?
data Conf = Conf {
        cf_index_list :: FilePath,
        cf_output     :: Maybe (BamMeta -> Iteratee [BamRec] IO ()),
        cf_stats_hdl  :: Handle,
        cf_num_stats  :: Int -> Int,
        cf_threshold  :: Double,
        cf_loudness   :: Loudness,
        cf_single     :: Bool,
        cf_samplesize :: Int,
        cf_readgroups :: [FilePath],
        cf_implied    :: [T.Text],
        cf_merge      :: Maybe ([U.Vector Nucleotides], [U.Vector Nucleotides]) }

defaultConf :: IO Conf
defaultConf = do ixdb <- getDataFileName "index_db.json"
                 return $ Conf {
                        cf_index_list = ixdb,
                        cf_output     = Nothing,
                        cf_stats_hdl  = stdout,
                        cf_num_stats  = \l -> max 20 $ l * 5 `div` 4,
                        cf_threshold  = 0.000005,
                        cf_loudness   = Normal,
                        cf_single     = False,
                        cf_samplesize = 50000,
                        cf_readgroups = [],
                        cf_implied    = [default_rgs],
                        cf_merge      = Nothing }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output"]          (ReqArg set_output   "FILE") "Send output to FILE",
    Option "I" ["index-database"]  (ReqArg set_index_db "FILE") "Read index database from FILE",
    Option "r" ["read-groups"]     (ReqArg set_rgs      "FILE") "Read read group definitions from FILE",
    Option "s" ["single-index"]    (NoArg           set_single) "Only consider first index",
    Option [ ] ["threshold"]       (ReqArg set_thresh   "FRAC") "Iterate till uncertainty is below FRAC",
    Option [ ] ["sample"]          (ReqArg set_sample    "NUM") "Sample NUM reads for mixture estimation",
    Option [ ] ["components"]      (ReqArg set_compo     "NUM") "Print NUM components of the mixture",
    Option [ ] ["nocontrol"]       (NoArg       set_no_control) "Suppress implied read groups for controls",
    Option "F" ["forward-adapter"] (ReqArg set_forward   "SEQ") "SEQ is a possible forward adapter",
    Option "R" ["reverse-adapter"] (ReqArg set_reverse   "SEQ") "SEQ is a possible reverse adapter",
    Option "v" ["verbose"]         (NoArg             set_loud) "Print more diagnostic messages",
    Option "q" ["quiet"]           (NoArg            set_quiet) "Print fewer diagnostic messages",
    Option "h?"["help", "usage"]   (NoArg        $ const usage) "Print this message and exit",
    Option "V" ["version"]         (NoArg        $  const vrsn) "Display version number and exit" ]
  where
    set_output  "-" c = return $ c { cf_output = Just $ pipeBamOutput, cf_stats_hdl = stderr }
    set_output   fp c = return $ c { cf_output = Just $ writeBamFile fp }
    set_index_db fp c = return $ c { cf_index_list = fp }
    set_rgs      fp c = return $ c { cf_readgroups = fp : cf_readgroups c }
    set_loud        c = return $ c { cf_loudness = Loud }
    set_quiet       c = return $ c { cf_loudness = Quiet }
    set_single      c = return $ c { cf_single = True }
    set_no_control  c = return $ c { cf_implied = [] }
    set_thresh    a c = readIO a >>= \x -> return $ c { cf_threshold = x }
    set_sample    a c = readIO a >>= \x -> return $ c { cf_samplesize = x }
    set_compo     a c = readIO a >>= \x -> return $ c { cf_num_stats = const x }

    set_forward   a c = readIO a >>= \x -> return $ c { cf_merge = Just $ first  (x:) $ fromMaybe ([],[]) $ cf_merge c }
    set_reverse   a c = readIO a >>= \x -> return $ c { cf_merge = Just $ second (x:) $ fromMaybe ([],[]) $ cf_merge c }

    usage = do pn <- getProgName
               putStrLn $ usageInfo ("Usage: " ++ pn ++ " [options] [bam-files]\n" ++
                                     "Decomposes a mix of libraries and assigns read groups.") options
               exitSuccess

    vrsn = do pn <- getProgName
              hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
              exitSuccess


adj_left :: Int -> Char -> L.Builder -> L.Builder
adj_left n c b = mconcat (replicate (n - fromIntegral (L.length t)) (L.singleton c)) <> L.fromLazyText t
  where t = L.toLazyText b

adj_left_text :: Int -> Char -> T.Text -> L.Builder
adj_left_text n c t = mconcat (replicate (n - T.length t) (L.singleton c)) <> L.fromText t

main :: IO ()
main = do
    (opts, files, errs) <- getOpt Permute options <$> getArgs
    unless (null errs) $ mapM_ (hPutStrLn stderr) errs >> exitFailure
    Conf{..} <- foldl (>>=) defaultConf opts
    when (null files) $ hPutStrLn stderr "no input files." >> exitFailure
    add_pg <- addPG $ Just version

    let notice  = case cf_loudness of Quiet -> \_ -> return () ; _ -> hPutStr stderr
        info    = case cf_loudness of Loud  -> hPutStr stderr ;  _ -> \_ -> return ()

    Both{..} <- B.readFile cf_index_list >>= \raw -> case decodeStrict' raw of
                    Nothing -> hPutStrLn stderr "Couldn't parse index database." >> exitFailure
                    Just  x | cf_single -> return $ x { p5is = single_placeholder }
                            | otherwise -> return   x

    rgdefs <- concatMap (readRGdefns (alias_names p7is) (alias_names p5is)) . (++) cf_implied <$> mapM T.readFile cf_readgroups
    notice $ "Got " ++ showNum (U.length (unique_indices p7is)) ++ " unique P7 indices and "
                    ++ showNum (U.length (unique_indices p5is)) ++ " unique P5 indices.\n"
    notice $ "Declared " ++ showNum (length rgdefs) ++ " read groups.\n"

    let n7     = U.length $ unique_indices p7is
        n5     = U.length $ unique_indices p5is
        stride = stride' n5
        vsize  = n7 * stride

    !rgs <- do let dup_error x y = error $ "Read groups " ++ show (fst x) ++ " and "
                                        ++ show (fst y) ++ " have the same indices."
               HM.fromListWith dup_error <$> sequence
                    [ VSM.replicate vsize (0::Double) >>= \dirt -> return ((i7,i5),(rg,dirt))
                    | RG !rg !i7 !i5 _ <- rgdefs ]

    let inspect = inspect' rgs (canonical_names p7is) (canonical_names p5is)

    ixvec <- concatInputs files >=> run $ gather cf_samplesize notice info
    notice $ "Got " ++ showNum (U.length ixvec) ++ " index pairs.\n"

    notice "decomposing mix "
    let go !n v = do v' <- iterEM ixvec (unique_indices p7is) (unique_indices p5is) v
                     case cf_loudness of Loud   -> hPutStrLn stderr [] >> inspect stderr 20 v'
                                         Normal -> hPutStr stderr "."
                                         Quiet  -> return ()
                     let d = VS.foldl' (\a -> max a . abs) 0 $ VS.zipWith (-) v v'
                     if n > 0 && d > cf_threshold * fromIntegral (U.length ixvec)
                          then go (n-1) v'
                          else do notice (if n == 0 then "\nmaximum number of iterations reached.\n"
                                                    else "\nmixture ratios converged.\n")
                                  return v'

    mix <- go (50::Int) $ naiveMix (U.length $ unique_indices p7is, U.length $ unique_indices p5is) (U.length ixvec)

    unlessQuiet cf_loudness $ do
            T.hPutStrLn cf_stats_hdl "\nfinal mixture estimate:"
            inspect cf_stats_hdl (cf_num_stats $ HM.size rgs) mix

    let maxlen = maximum $ map (B.length . rgid) rgdefs
        ns7 = canonical_names p7is
        ns5 = canonical_names p5is
        num = 7

    case cf_output of
        Nothing  -> do  unlessQuiet cf_loudness $ do
                            T.hPutStrLn cf_stats_hdl "\nmaximum achievable quality, top pollutants:"
                            forM_ (sortOn (fst.snd) $ HM.toList rgs) $ \((i7,i5), (rgid,_)) -> do
                                (p,_,_) <- class1 HM.empty (unique_indices p7is) (unique_indices p5is) mix
                                                  (unique_indices p7is U.! i7, unique_indices p5is U.! i5)

                                let qmax = negate . round $ 10 / log 10 * log p :: Int
                                L.hPutStrLn cf_stats_hdl . L.toLazyText $
                                        adj_left_text maxlen ' ' (T.decodeUtf8 rgid) <>
                                        L.fromText ": " <>
                                        adj_left 4 ' ' (L.singleton 'Q' <> L.decimal (max 0 qmax))

        Just out -> do  concatInputs files >=> run $ \hdr ->
                            let hdr' = hdr { meta_other_shit =
                                              [ os | os@(k,_) <- meta_other_shit hdr, k /= "RG" ] ++
                                              HM.elems (HM.fromList [ (rgid, ("RG", ("ID",rgid):tags)) | RG{..} <- rgdefs ] ) }
                            in mapStreamM (\br -> do
                                    let b = unpackBam br
                                        eff_rgs | not (isPaired b) = rgs
                                                | isFirstMate b    = rgs
                                                | otherwise        = HM.empty
                                    (p,i7,i5) <- class1 eff_rgs (unique_indices p7is) (unique_indices p5is) mix
                                                                (fromTags "XI" "YI" br, fromTags "XJ" "YJ" br)
                                    let q = negate . round $ 10 / log 10 * log p
                                        ex = deleteE "ZR" . deleteE "Z0" . deleteE "Z2" . updateE "Z1" (Int q) $
                                             updateE "ZX" (Text $ T.encodeUtf8 $ T.concat [ ns7 V.! i7, ",", ns5 V.! i5 ]) $
                                             case HM.lookup (i7,i5) rgs of
                                               Nothing      -> deleteE "RG" $ b_exts b
                                               Just (rgn,_) -> updateE "RG" (Text rgn) $ b_exts b
                                    return $ case lookup "ZQ" ex of
                                                Just (Text t) | BS.null t' -> b { b_exts = deleteE "ZQ" ex
                                                                                , b_flag = b_flag b .&. complement flagFailsQC }
                                                              | otherwise  -> b { b_exts = updateE "ZQ" (Text t') ex }
                                                  where
                                                    t' = BS.filter (\c -> c /= 'C' && c /= 'I' && c /= 'W') t
                                                _                          -> b { b_exts = ex
                                                                                , b_flag = b_flag b .&. complement flagFailsQC }) =$
                               progressNum "writing " 0x100000 info =$
                               maybe (mergeTrimBam default_fwd_adapters default_rev_adapters)
                                     (uncurry mergeTrimBam) cf_merge =$
                               out (add_pg hdr')

                        unlessQuiet cf_loudness $ do
                            grand_total <- foldM (\ !acc (_,dirt) -> VS.freeze dirt >>= return . (+) acc . VS.sum) 0 (HM.elems rgs)
                            T.hPutStrLn cf_stats_hdl "\nmaximum achievable and average quality, top pollutants:"
                            forM_ (sortOn (fst.snd) $ HM.toList rgs) $ \((i7,i5), (rgid,dirt_)) -> do
                                dirt <- VS.freeze dirt_
                                (p,_,_) <- class1 HM.empty (unique_indices p7is) (unique_indices p5is) mix
                                                  (unique_indices p7is U.! i7, unique_indices p5is U.! i5)

                                let total  = VS.sum dirt
                                    others = VS.sum $ VS.ifilter (\i _ -> i /= i7 * stride + i5) dirt
                                    qmax = negate . round $ 10 / log 10 * log p :: Int
                                    qavg = negate . round $ 10 / log 10 * log (others/total) :: Int

                                v <- U.unsafeThaw . U.fromListN (VS.length dirt) . zip [0..] . VS.toList $ dirt
                                V.sortBy (\(_,a) (_,b) -> compare b a) v -- meh.
                                v' <- U.unsafeFreeze v

                                let fmt_one (i,n) =
                                        let (i7', i5') = i `quotRem` stride
                                            chunk = L.formatRealFloat L.Fixed (Just 2) (100*n/total) <> L.singleton '%' <>
                                                    L.singleton ' ' <> L.fromText (ns7 V.! i7') <>
                                                    L.singleton '/' <> L.fromText (ns5 V.! i5') <>
                                                    case HM.lookup (i7',i5') rgs of
                                                        Nothing     -> mempty
                                                        Just (rg,_) -> L.singleton ' ' <> L.singleton '(' <>
                                                                       L.fromText (T.decodeUtf8 rg) <> L.singleton ')'
                                        in if (i7 == i7' && i5 == i5') || i5' >= n5 then id else (:) chunk

                                when (total >= 1) . L.hPutStrLn cf_stats_hdl . L.toLazyText $
                                        adj_left_text maxlen ' ' (T.decodeUtf8 rgid) <>
                                        L.singleton ':' <> L.singleton ' ' <>
                                        adj_left 4 ' ' (L.singleton 'Q' <> L.decimal (max 0 qmax)) <> L.fromText ", " <>
                                        adj_left 4 ' ' (L.singleton 'Q' <> L.decimal (max 0 qavg)) <> L.fromText ", " <>
                                        L.fromString (showNum (round total :: Int)) <> L.fromText " (" <>
                                        L.formatRealFloat L.Fixed (Just 2) (100*total/grand_total) <> L.fromText "%); " <>
                                        foldr1 (\a b -> a <> L.fromText ", " <> b)
                                            (take num $ U.foldr fmt_one [] v')

inspect' :: HM.HashMap (Int,Int) (B.ByteString, t) -> V.Vector T.Text -> V.Vector T.Text -> Handle -> Int -> Mix -> IO ()
inspect' rgs n7 n5 hdl num_ mix = do
    -- Due to padding, we get invalid indices here.  Better filter them
    -- out, because we sure can't print them later.
    let num = min num_ $ V.length n5 * V.length n7
    v <- U.unsafeThaw $ U.fromListN (V.length n5 * V.length n7) $
                filter (\(i,_) -> i `rem` stride' (V.length n5) < V.length n5) $
                zip [0..] $ VS.toList mix
    V.partialSortBy (\(_,a) (_,b) -> compare b a) v num         -- meh.
    v' <- U.unsafeFreeze v

    let total  = U.sum . U.map snd $ v'
        others = U.sum . U.map snd . U.drop num $ v'

    U.forM_ (U.take num v') $ \(i,n) -> do
       let (i7, i5) = i `quotRem` stride' (V.length n5)
       L.hPutStrLn hdl . L.toLazyText $
            adj_left_text 7 ' ' (n7 V.! i7) <> L.singleton ',' <> L.singleton ' ' <>
            adj_left_text 7 ' ' (n5 V.! i5) <> L.singleton ':' <> L.singleton ' ' <>
            adj_left 8 ' ' (L.formatRealFloat L.Fixed (Just 3) (100 * n / total)) <> L.singleton '%' <> L.singleton ' ' <>
            case HM.lookup (i7,i5) rgs of
                Nothing     -> mempty
                Just (rg,_) -> L.singleton '(' <> L.fromText (T.decodeUtf8 rg) <> L.singleton ')'

    L.hPutStrLn hdl . L.toLazyText $
        L.fromLazyText "      all others: " <>
        adj_left 8 ' ' (L.formatRealFloat L.Fixed (Just 3) (100 * others / total)) <>
        L.singleton '%'

