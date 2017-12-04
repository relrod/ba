{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE PatternSynonyms #-}
-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types
-- Copyright : (C) 2017 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- We emphasize making illegal states unrepresentable. As such, creating
-- something representing a particular kind of
-- <http://rosalind.info/glossary/genetic-string/ genetic string> will be done
-- by way of smart constructors. The default constructors will not be exported.
--
-- We represent three different kinds of genetic strings (as per Rosalind),
-- 'DNA', 'RNA', and 'Protein'.
--
-- A 'DNA' is created with 'mkDNA', an 'RNA' is created with 'mkRNA', and a
-- 'Protein' is created with 'mkProtein'.
--
-- Because we don\'t export the constructors, pattern matching becomes
-- impossible. To rectify this, we use PatternSynonyms and export patterns for
-- each type. These are 'PDNA', 'PRNA', and 'PProtein', respectively.
-- Internally, we represent sequences as lazy 'T.Text's, and this is what you
-- deconstruct one using one of these pattern synonyms.
----------------------------------------------------------------------------

module Bio.Algorithm.Types
  (
  -- * Macromolecule types
    DNA
  , RNA
  -- * Prisms (DNA)
  , stringDNA
  , lazyTextDNA
  , strictTextDNA
  , lazyBSDNA
  , strictBSDNA
  -- * Prisms (RNA)
  , stringRNA
  , lazyTextRNA
  , strictTextRNA
  , lazyBSRNA
  , strictBSRNA
  -- * Transcription isomorphism
  , transcribe
  , reverseTranscribe
  -- * Trivial isomorphisms
  , dnaIso
  , rnaIso
  -- * Helper Functions
  , mapDNA
  , mapRNA
  , isDNAChar
  , isRNAChar
  ) where

import Control.Lens
import qualified Data.ByteString.Char8 as C8
import qualified Data.ByteString.Lazy.Char8 as LC8
import Data.Semigroup
import qualified Data.Text as T
import qualified Data.Text.Encoding as TE
import qualified Data.Text.Lazy as TL
import qualified Data.Text.Lazy.Encoding as TLE

-- | A 'DNA' sequence. This is constructed using the 'mkDNA' smart constructor
-- below.
newtype DNA = DNA TL.Text deriving (Eq, Ord, Show, Monoid, Semigroup)

-- | An 'RNA' sequence. This is constructed using the 'mkRNA' smart constructor
-- below.
newtype RNA = RNA TL.Text deriving (Eq, Ord, Show, Monoid, Semigroup)

-- | There a many different kinds of structures which can be converted into a
-- 'DNA'. We specify these using 'Prism''s.
class AsDNA a where
  _DNA :: Prism' a DNA

-- | Map over a 'DNA' sequence if possible. In general, prefer to use more
-- \"Lensy\" ways of doing this.
mapDNA :: (AsDNA s, AsDNA t) => (s -> t) -> DNA -> Maybe DNA
mapDNA f dna = dna ^. re _DNA . to f ^? _DNA

-- | Map over an 'RNA' sequence if possible. In general, prefer to use more
-- \"Lensy\" ways of doing this.
mapRNA :: (AsRNA s, AsRNA t) => (s -> t) -> RNA -> Maybe RNA
mapRNA f rna = rna ^. re _RNA . to f ^? _RNA

-- | Helper function for helping to validate a valid 'DNA' string.
isDNAChar :: Char -> Bool
isDNAChar 'A' = True
isDNAChar 'C' = True
isDNAChar 'G' = True
isDNAChar 'T' = True
isDNAChar _   = False

-- | DNA is, of course, isomorphic to itself under the identity function.
dnaIso :: Iso' DNA DNA
dnaIso = iso id id

-- | Free since every 'Iso' is a 'Prism'.
instance AsDNA DNA where
  _DNA = dnaIso

stringDNA :: Prism' String DNA
stringDNA = prism' (\(DNA tl) -> TL.unpack tl) toDNA
  where
    toDNA s =
      if length (filter isDNAChar s) == length s
      then Just (DNA (TL.pack s))
      else Nothing

instance AsDNA String where
  _DNA = stringDNA

lazyTextDNA :: Prism' TL.Text DNA
lazyTextDNA = prism' (\(DNA tl) -> tl) toDNA
  where
    toDNA s =
      if TL.length (TL.filter isDNAChar s) == TL.length s
      then Just (DNA s)
      else Nothing

instance AsDNA TL.Text where
  _DNA = lazyTextDNA

strictTextDNA :: Prism' T.Text DNA
strictTextDNA = prism' (\(DNA t) -> TL.toStrict t) (\t -> t ^? lazy . _DNA)

instance AsDNA T.Text where
  _DNA = strictTextDNA

lazyBSDNA :: Prism' LC8.ByteString DNA
lazyBSDNA = prism' (\(DNA tl) -> TLE.encodeUtf8 tl) toDNA
  where
    toDNA s =
      if LC8.length (LC8.filter isDNAChar s) == LC8.length s
      then Just (DNA (TLE.decodeUtf8 s))
      else Nothing

instance AsDNA LC8.ByteString where
  _DNA = lazyBSDNA

strictBSDNA :: Prism' C8.ByteString DNA
strictBSDNA =
  prism' (\(DNA t) -> t ^. strict . to TE.encodeUtf8) (\t -> t ^? lazy . _DNA)

instance AsDNA C8.ByteString where
  _DNA = strictBSDNA

-- | There a many different kinds of structures which can be converted into a
-- 'RNA'. We specify these using 'Prism''s.
class AsRNA a where
  _RNA :: Prism' a RNA

-- | Helper function for helping to validate a valid 'RNA' string.
isRNAChar :: Char -> Bool
isRNAChar 'A' = True
isRNAChar 'C' = True
isRNAChar 'G' = True
isRNAChar 'U' = True
isRNAChar _   = False

-- | RNA is, of course, isomorphic to itself under the identity function.
rnaIso :: Iso' RNA RNA
rnaIso = iso id id

-- | Free since every 'Iso' is a 'Prism'.
instance AsRNA RNA where
  _RNA = rnaIso

stringRNA :: Prism' String RNA
stringRNA = prism' (\(RNA tl) -> TL.unpack tl) toRNA
  where
    toRNA s =
      if length (filter isRNAChar s) == length s
      then Just (RNA (TL.pack s))
      else Nothing

instance AsRNA String where
  _RNA = stringRNA

lazyTextRNA :: Prism' TL.Text RNA
lazyTextRNA = prism' (\(RNA tl) -> tl) toRNA
  where
    toRNA s =
      if TL.length (TL.filter isRNAChar s) == TL.length s
      then Just (RNA s)
      else Nothing

instance AsRNA TL.Text where
  _RNA = lazyTextRNA

strictTextRNA :: Prism' T.Text RNA
strictTextRNA = prism' (\(RNA t) -> TL.toStrict t) (\t -> t ^? lazy . _RNA)

instance AsRNA T.Text where
  _RNA = strictTextRNA

lazyBSRNA :: Prism' LC8.ByteString RNA
lazyBSRNA = prism' (\(RNA tl) -> TLE.encodeUtf8 tl) toRNA
  where
    toRNA s =
      if LC8.length (LC8.filter isRNAChar s) == LC8.length s
      then Just (RNA (TLE.decodeUtf8 s))
      else Nothing

instance AsRNA LC8.ByteString where
  _RNA = lazyBSRNA

strictBSRNA :: Prism' C8.ByteString RNA
strictBSRNA =
  prism' (\(RNA t) -> t ^. strict . to TE.encodeUtf8) (\t -> t ^? lazy . _RNA)

instance AsRNA C8.ByteString where
  _RNA = strictBSRNA

-- | There is an isomorphism between 'DNA' and 'RNA' defined by 'RNA'
-- transcription and reverse-transcription.
transcribe :: Iso' DNA RNA
transcribe = iso toRNA fromRNA
  where
    comp 'A' = 'U'
    comp 'C' = 'G'
    comp 'G' = 'C'
    comp 'T' = 'A'
    comp x = x -- should never happen
    revComp 'U' = 'A'
    revComp 'G' = 'C'
    revComp 'C' = 'G'
    revComp 'A' = 'T'
    revComp x = x -- should never happen
    toRNA (DNA dna) = RNA (TL.map comp dna)
    fromRNA (RNA rna) = DNA (TL.map revComp rna)

-- | A helper to reverse-transcribe back to 'DNA'.
--
-- @reverseTranscribe = from transcribe@
reverseTranscribe :: Iso' RNA DNA
reverseTranscribe = from transcribe
