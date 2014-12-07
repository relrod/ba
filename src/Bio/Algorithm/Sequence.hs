module Bio.Algorithm.Sequence where

import Control.Arrow
import Data.List
import qualified Data.Map as Map

kmers :: Ord a => Int -> [a] -> [[a]]
kmers k s = getMaxes . frequency . kmers'  k $ s
  where
    frequency :: Ord a => [a] -> Map.Map a Int
    frequency = Map.fromList . fmap (head &&& length) . group . sort

    getMaxes :: Ord a => Map.Map k a -> [k]
    getMaxes m = Map.keys (Map.filter (==f) m)
      where
        f = maximum . Map.elems $ m

    kmers' :: Int -> [a] -> [[a]]
    kmers' _ [] = []
    kmers' k' s' = if length s' >= k'
                   then take k' s' : kmers' k' (drop 1 s')
                   else kmers' k' (drop 1 s')

reverseComplement :: String -> String
reverseComplement = reverse . map complement
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
