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
, ltClumps
) where

import Bio.Algorithm.Types
import Control.Applicative
import Control.Arrow
import Control.Monad
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

-- | Find and list all distinct (L, t)-clumps of the given length in the given
-- string.
--
-- This is done by scanning the entire genome with a window of size @l@. In each
-- window, we get all @k@-mers, see if they exist @t@ times, and if they do, we
-- add them to the result.
ltClumps :: Int           -- ^ @k@. The length of the k-mers (i.e., the /k/ of the "k-mer")
         -> Int           -- ^ @l@. The interval size to scan
         -> Int           -- ^ @t@. How many times it must appear
         -> RawSequence   -- ^ The sequence to search in
         -> [RawSequence] -- ^ The resulting list of k-mers
ltClumps k l t (RawSequence s) = nub . tTimes . map (kmers k . RawSequence) $ windows s
  where
    windows "" = []
    windows s' = BL.take (fromIntegral l) s' : windows (BL.drop 1 s')

    frequency :: Ord a => [a] -> [(a, Int)]
    frequency = fmap (head &&& length) . group . sort

    tTimes x = map fst . filter (\(_, y) -> y == t) . join . map frequency $ x

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
