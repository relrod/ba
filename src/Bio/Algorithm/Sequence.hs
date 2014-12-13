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

-- $setup
-- >>> import Bio.Algorithm.Types
-- >>> import qualified Data.ByteString.Lazy.Char8 as BL
-- >>> import qualified Data.Text as T
-- >>> import Control.Lens

-- | Find and list all of the kmers of the given length in the given string.
--
-- >>> kmers 3 "AGATCGAGTG"
-- ["AGA","AGT","ATC","CGA","GAG","GAT","GTG","TCG"]
--
-- >>> kmers 5 "GTAGAGCTGT"
-- ["AGAGC","AGCTG","GAGCT","GCTGT","GTAGA","TAGAG"]
kmers :: Ord a =>
         Int   -- ^ The length of the k-mers (i.e., the /k/ of the "k-mer")
      -> [a]   -- ^ The string (or list) to search in
      -> [[a]] -- ^ The resulting list of k-mers
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

-- | Return the reverse complement for some DNA sequence.
--
-- >>> dnaReverseComplement (RawSequence (BL.pack "CCTTGGAA"))
-- RawSequence "TTCCAAGG"
--
-- >>> T.pack "CCTTGGAA" ^. lazy . _RawSequence . to dnaReverseComplement
-- RawSequence "TTCCAAGG"
dnaReverseComplement :: RawSequence -> RawSequence
dnaReverseComplement (RawSequence s) = RawSequence . BL.reverse . BL.map complement $ s
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

-- | Return the reverse complement for some RNA sequence.
--
-- >>> rnaReverseComplement (RawSequence (BL.pack "CCTTGGAA"))
-- RawSequence "UUCCAAGG"
--
-- >>> T.pack "CCTTGGAA" ^. lazy . _RawSequence . to rnaReverseComplement
-- RawSequence "UUCCAAGG"
rnaReverseComplement :: RawSequence -> RawSequence
rnaReverseComplement (RawSequence s) = RawSequence . BL.reverse . BL.map complement $ s
  where
    complement 'A' = 'U'
    complement 'a' = 'u'
    complement 'U' = 'A'
    complement 'u' = 'a'
    complement 'C' = 'G'
    complement 'c' = 'g'
    complement 'G' = 'C'
    complement 'g' = 'c'
    complement x   = x
