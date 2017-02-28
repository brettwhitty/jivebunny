{-# LANGUAGE GeneralizedNewtypeDeriving, TemplateHaskell #-}
{-# LANGUAGE TypeFamilies, MultiParamTypeClasses, CPP    #-}
module Index where

-- ^ This tiny module defines the 'Index' type and derives the 'Unbox'
-- instance.  That dramatically lowers the chance that template haskell
-- runs into problems :(

import Bio.Prelude
import Data.Vector.Unboxed.Deriving

#if __GLASGOW_HASKELL__ == 704
import Data.Vector.Generic          ( Vector(..) )
import Data.Vector.Generic.Mutable  ( MVector(..) )
#endif

-- | An index sequence must have at most eight bases.  We represent a
-- base and its quality score in a single byte:  the top three bits are
-- the base ("ACGTN" = [0,1,3,2,7]), the lower five bits are the quality,
-- clamped to 31.

newtype Index = Index Word64 deriving (Storable, Eq)

instance Hashable Index where
    hashWithSalt salt (Index x) = hashWithSalt salt x
    hash (Index x) = hash x

instance Show Index where
    show (Index x) = [ "ACTGNNNN" !! fromIntegral b | i <- [56,48..0], let b = (x `shiftR` (i+5)) .&. 0x7 ]
            ++ 'q' : [ chr (fromIntegral q+33)      | i <- [56,48..0], let q = (x `shiftR` i) .&. 0x1F ]

derivingUnbox "Index" [t| Index -> Word64 |] [| \ (Index i) -> i |] [| Index |]


