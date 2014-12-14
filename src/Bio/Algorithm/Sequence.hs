{-# LANGUAGE OverloadedStrings #-}
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
module Bio.Algorithm.Sequence (
  -- * Complements and Reverse Complements
  dnaReverseComplement
, rnaReverseComplement

  -- * Transcription
, dnaToRna

  -- * k-mer algorithms
, kmers
, maxKmers
) where

import Bio.Algorithm.Types
import Control.Applicative
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
-- >>> kmers 3 (RawSequence $ BL.pack "AGATCGAGTG")
-- [RawSequence "AGA",RawSequence "GAT",RawSequence "ATC",RawSequence "TCG",RawSequence "CGA",RawSequence "GAG",RawSequence "AGT",RawSequence "GTG"]
--
-- >>> kmers 5 (RawSequence $ BL.pack "GTAGAGCTGT")
-- [RawSequence "GTAGA",RawSequence "TAGAG",RawSequence "AGAGC",RawSequence "GAGCT",RawSequence "AGCTG",RawSequence "GCTGT"]
--
-- >>> "gatatat" ^. _RawSequence . to (kmers 4)
-- [RawSequence "gata",RawSequence "atat",RawSequence "tata",RawSequence "atat"]
kmers :: Int           -- ^ The length of the k-mers (i.e., the /k/ of the "k-mer")
      -> RawSequence   -- ^ The sequence to search in
      -> [RawSequence] -- ^ The resulting list of k-mers
kmers k (RawSequence s) = RawSequence <$> kmers' k s
  where
    kmers' :: Int -> BL.ByteString -> [BL.ByteString]
    kmers' _ "" = []
    kmers' k' s' = if BL.length s' >= fromIntegral k'
                   then BL.take (fromIntegral k') s' : kmers' k' (BL.drop 1 s')
                   else kmers' k' (BL.drop 1 s')

-- | Return only the kmers that appear the highest number of times.
--
-- >>> "gatatat" ^. _RawSequence . to (maxKmers . kmers 4)
-- [RawSequence "atat"]
maxKmers :: [RawSequence] -> [RawSequence]
maxKmers = getMaxes . frequency
  where
    frequency :: Ord a => [a] -> Map.Map a Int
    frequency = Map.fromList . fmap (head &&& length) . group . sort

    getMaxes :: (Ord a, Ord k) => Map.Map k a -> [k]
    getMaxes m = Map.keys (Map.filter (==f) m)
      where
        f = maximum . Map.elems $ m

-- | Return the reverse complement for some DNA sequence.
--
-- >>> dnaReverseComplement (DNA (RawSequence (BL.pack "CCTTGGAA")))
-- DNA (RawSequence "TTCCAAGG")
--
-- >>> T.pack "CCTTGGAA" ^. lazy . _RawSequence . to (dnaReverseComplement . DNA)
-- DNA (RawSequence "TTCCAAGG")
dnaReverseComplement :: DNA -> DNA
dnaReverseComplement (DNA (RawSequence s)) = DNA . RawSequence . BL.reverse . BL.map complement $ s
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
-- >>> rnaReverseComplement (RNA (RawSequence (BL.pack "GAUGGAACUUGACUACGUAAAUU")))
-- RNA (RawSequence "AAUUUACGUAGUCAAGUUCCAUC")
--
-- >>> T.pack "GAUGGAACUUGACUACGUAAAUU" ^. lazy . _RawSequence . to (rnaReverseComplement . RNA)
-- RNA (RawSequence "AAUUUACGUAGUCAAGUUCCAUC")
rnaReverseComplement :: RNA -> RNA
rnaReverseComplement (RNA (RawSequence s)) = RNA . RawSequence . BL.reverse . BL.map complement $ s
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

-- | Transcribes a DNA sequence into RNA.
--
-- Transcribe a DNA sequence:
--
-- >>> BL.pack "GATGGAACTTGACTACGTAAATT" ^. _RawSequence . to (dnaToRna . DNA)
-- RNA (RawSequence "GAUGGAACUUGACUACGUAAAUU")
--
-- Transcribe and then take the 'rnaReverseComplement':
--
-- >>> BL.pack "GATGGAACTTGACTACGTAAATT" ^. _RawSequence . to (rnaReverseComplement . dnaToRna . DNA)
-- RNA (RawSequence "AAUUUACGUAGUCAAGUUCCAUC")
dnaToRna :: DNA -> RNA
dnaToRna (DNA (RawSequence s)) = RNA . RawSequence . BL.map rna $ s
  where
    rna 'T' = 'U'
    rna 't' = 'u'
    rna x   = x
