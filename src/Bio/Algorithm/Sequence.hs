-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Sequence
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- Algorithms dealing with bioinformatics sequences. Right now it doesn't do
-- much.
----------------------------------------------------------------------------
module Bio.Algorithm.Sequence where

import Bio.Algorithm.Types
import Control.Arrow
import qualified Data.ByteString.Lazy.Char8 as BL
import Data.List
import qualified Data.Map as Map

kmers :: Ord a => Int -> [a] -> [[a]]
kmers k s = getMaxes . frequency . kmers'  k $ s
  where
    frequency :: Ord a => [a] -> Map.Map a Int
    frequency = Map.fromList . fmap (head &&& length) . group . sort

    getMaxes :: (Ord a, Ord k) => Map.Map k a -> [k]
    getMaxes m = Map.keys (Map.filter (==f) m)
      where
        f = maximum . Map.elems $ m

    kmers' :: Int -> [a] -> [[a]]
    kmers' _ [] = []
    kmers' k' s' = if length s' >= k'
                   then take k' s' : kmers' k' (drop 1 s')
                   else kmers' k' (drop 1 s')

reverseComplement :: RawSequence -> RawSequence
reverseComplement (RawSequence s) = RawSequence . BL.reverse . BL.map complement $ s
  where
    complement 'A' = 'T'
    complement 'a' = 't'
    complement 'T' = 'A'
    complement 't' = 'a'
    complement 'C' = 'G'
    complement 'c' = 'g'
    complement 'G' = 'C'
    complement 'g' = 'c'
    complement x   = x
