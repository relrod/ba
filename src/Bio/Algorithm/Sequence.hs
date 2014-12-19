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

  -- * k-mer algorithms
, kmers
, maxKmers
, ltClumps
, approximateMatchIndices

  -- * RNA -> Protein Translation
, translate

  -- * Skew
, skews
) where

import Bio.Algorithm.Types
import Control.Applicative
import Control.Arrow
import Control.Lens
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
--
-- >>> "gatatat" ^. _RawSequence . to DNA . to (kmers 4)
-- [DNA (RawSequence "gata"),DNA (RawSequence "atat"),DNA (RawSequence "tata"),DNA (RawSequence "atat")]
--
-- >>> "gatatat" ^. _RawSequence . to RNA . to (kmers 4)
-- [RNA (RawSequence "gata"),RNA (RawSequence "atat"),RNA (RawSequence "tata"),RNA (RawSequence "atat")]
kmers :: AsRawSequence s
      => Int -- ^ The length of the k-mers (i.e., the /k/ of the "k-mer")
      -> s   -- ^ The sequence to search in
      -> [s] -- ^ The resulting list of k-mers
kmers k s = (_RawSequence #) . RawSequence <$> kmers' k (_RawSequence # (s ^. _RawSequence))
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
--
-- >>> "gatatat" ^. _RawSequence . to (maxKmers . kmers 4 . DNA)
-- [DNA (RawSequence "atat")]
--
-- >>> "gatatat" ^. _RawSequence . to (maxKmers . kmers 4 . RNA)
-- [RNA (RawSequence "atat")]
maxKmers :: (AsRawSequence s, Ord s)
         => [s]
         -> [s]
maxKmers s = getMaxes . frequency $ fmap ((_RawSequence #) . (^. _RawSequence)) s
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
--
-- >>> "gatcagcataagggtccctgcaatgcatgacaagcctgcagttgttttac" ^. _RawSequence . to (ltClumps 4 25 3)
-- [RawSequence "tgca"]
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

-- | Find skews (difference in number of G and number of C occurances in some
-- sequence) by removing a nucleotide from the right side and recursing over the
-- remaining sequence.
--
-- Can be used to solve the Rosalind "Minimum Skew Problem" like this:
--  @minSkewIndices s = map (+1) . findIndices (== minimum (skews s)) $ (skews s)@
-- although it is too slow to actually be accepted.
skews :: AsRawSequence s => s -> [Int]
skews s = reverse . map fromIntegral $ runSkews (_RawSequence # (s ^. _RawSequence))
  where
    runSkews "" = []
    runSkews s' =
      let
        numG = BL.length . BL.filter (liftM2 (||) (== 'G') (== 'g')) $ s'
        numC = BL.length . BL.filter (liftM2 (||) (== 'C') (== 'c')) $ s'
      in (numG - numC) : runSkews (BL.init s')

-- | Finds approximate matches for some given string, in a given string.
--
-- We define an approximate match as a match where the hamming distance is less
-- than the given @d@.
--
-- >>> approximateMatchIndices (BL.pack "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC") (BL.pack "ATTCTGGA") 3
-- [6,7,26,27,78]
approximateMatchIndices :: BL.ByteString -- ^ The string to search in.
                        -> BL.ByteString -- ^ The string to search for.
                        -> Int           -- ^ @d@. Max number of mutations.
                        -> [Int]         -- ^ Indices of such matches.
approximateMatchIndices i pat d =
  runApproximation $ filter ((BL.length pat ==) . BL.length) . windows $ i
  where
    hamming x y = length . filter (uncurry (/=)) $ zip x y

    windows "" = []
    windows s' = BL.take (BL.length pat) s' : windows (BL.drop 1 s')

    runApproximation = findIndices (\x -> hamming (BL.unpack x) (BL.unpack pat) <= d)

-- | Translate t'RNA' to 'Protein'.
--
-- Returns 'Nothing' if given an invalid RNA sequence.
-- Stops immediately and returns the result so far, upon encountering a stop
-- codon ("UAA", "UAG", "UGA").
--
-- >>> translate (RNA ("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGAAUG" ^. _RawSequence))
-- Just [Protein 'M',Protein 'A',Protein 'M',Protein 'A',Protein 'P',Protein 'R',Protein 'T',Protein 'E',Protein 'I',Protein 'N',Protein 'S',Protein 'T',Protein 'R',Protein 'I',Protein 'N',Protein 'G']
--
-- >>> translate (RNA ("AUGGCfoobarfoobarCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGAAUG" ^. _RawSequence))
-- Nothing
--
-- >>> translate (RNA ("GAAUAAfoobarfoobar" ^. _RawSequence))
-- Just [Protein 'E']
translate :: RNA -> Maybe [Protein]
translate (RNA (RawSequence s)) = fmap (fmap Protein) (translate' s)
  where
    translate' "" = Just []
    translate' b = do
      p <- rnaToProtein (BL.take 3 b)
      case p of
       Stop -> return []
       Protein protein -> do
         x <- translate' (BL.drop 3 b)
         return (protein : x)

    rnaToProtein b = Map.lookup b proteinMap
      where
        a's   = ["GCU", "GCC", "GCA", "GCG"]
        c's   = ["UGU", "UGC"]
        d's   = ["GAU", "GAC"]
        e's   = ["GAA", "GAG"]
        f's   = ["UUU", "UUC"]
        g's   = ["GGU", "GGC", "GGA", "GGG"]
        h's   = ["CAU", "CAC"]
        i's   = ["AUU", "AUC", "AUA"]
        k's   = ["AAA", "AAG"]
        l's   = ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"]
        m's   = ["AUG"]
        n's   = ["AAU", "AAC"]
        p's   = ["CCU", "CCC", "CCA", "CCG"]
        q's   = ["CAA", "CAG"]
        r's   = ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"]
        s's   = ["AGU", "AGC", "UCU", "UCC", "UCA", "UCG"]
        t's   = ["ACU", "ACC", "ACA", "ACG"]
        v's   = ["GUU", "GUC", "GUA", "GUG"]
        w's   = ["UGG"]
        y's   = ["UAU", "UAC"]
        stops = ["UAA", "UAG", "UGA"]

        proteinMap = Map.fromList $ join [
            map (flip (,) (Protein 'A')) a's
          , map (flip (,) (Protein 'C')) c's
          , map (flip (,) (Protein 'D')) d's
          , map (flip (,) (Protein 'E')) e's
          , map (flip (,) (Protein 'F')) f's
          , map (flip (,) (Protein 'G')) g's
          , map (flip (,) (Protein 'H')) h's
          , map (flip (,) (Protein 'I')) i's
          , map (flip (,) (Protein 'K')) k's
          , map (flip (,) (Protein 'L')) l's
          , map (flip (,) (Protein 'M')) m's
          , map (flip (,) (Protein 'N')) n's
          , map (flip (,) (Protein 'P')) p's
          , map (flip (,) (Protein 'Q')) q's
          , map (flip (,) (Protein 'R')) r's
          , map (flip (,) (Protein 'S')) s's
          , map (flip (,) (Protein 'T')) t's
          , map (flip (,) (Protein 'V')) v's
          , map (flip (,) (Protein 'W')) w's
          , map (flip (,) (Protein 'Y')) y's
          , map (flip (,) Stop) stops
          ]
