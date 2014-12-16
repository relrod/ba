{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE OverloadedStrings #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types.Sequence
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- RNA sequences
----------------------------------------------------------------------------
module Bio.Algorithm.Types.Sequence where

import Bio.Algorithm.Types.RawSequence
import Control.Lens
import qualified Data.ByteString.Lazy.Char8 as BL

-- $setup
-- >>> import Bio.Algorithm.Sequence

----------------------------------------------------------------------------
-- DNA
----------------------------------------------------------------------------
newtype DNA = DNA RawSequence deriving (Eq, Ord, Show)

instance AsRawSequence DNA where
  _RawSequence = iso (\(DNA r) -> r) DNA
  {-# INLINE _RawSequence #-}

----------------------------------------------------------------------------
-- RNA
----------------------------------------------------------------------------
newtype RNA = RNA RawSequence deriving (Eq, Ord, Show)

instance AsRawSequence RNA where
  _RawSequence = iso (\(RNA r) -> r) RNA
  {-# INLINE _RawSequence #-}

-- | Transcribes a DNA sequence into RNA.
--
-- There exist (several, but one useful) isomorphism between 'DNA' and 'DNA' as
-- defined by David Spivak\'s
-- <http://math.mit.edu/~dspivak/teaching/sp13/CT4S--static.pdf Category Theory for Scientists>.
--
-- We provide that isomorphism here.
--
-- Transcribe a DNA sequence:
--
-- >>> BL.pack "GATGGAACTTGACTACGTAAATT" ^. _RawSequence . to DNA . from rnaDnaIso
-- RNA (RawSequence "GAUGGAACUUGACUACGUAAAUU")
--
-- >>> BL.pack "GAUGGAACUUGACUACGUAAAUU" ^. _RawSequence . to RNA . rnaDnaIso
-- DNA (RawSequence "GATGGAACTTGACTACGTAAATT")
--
-- Transcribe and then take the 'rnaReverseComplement':
--
-- λ> BL.pack "GATGGAACTTGACTACGTAAATT" ^. _RawSequence . to DNA . from rnaDnaIso . to rnaReverseComplement
-- RNA (RawSequence "AAUUUACGUAGUCAAGUUCCAUC")
rnaDnaIso :: Iso' RNA DNA
rnaDnaIso = iso transcribeRnaToDna transcribeDnaToRna
  where
    transcribeRnaToDna (RNA (RawSequence s)) = (BL.map rnaToDna s) ^. _RawSequence . to DNA
    transcribeDnaToRna (DNA (RawSequence s)) = (BL.map dnaToRna s) ^. _RawSequence . to RNA

    rnaToDna 'U' = 'T'
    rnaToDna 'u' = 't'
    rnaToDna x   = x

    dnaToRna 'T' = 'U'
    dnaToRna 't' = 'u'
    dnaToRna x   = x

-- | Helper to make using 'rnaDnaIso' when converting from 'DNA' to 'RNA'.
--
-- >>> rna # DNA ("GATGGAACTTGACTACGTAAATT" ^. _RawSequence)
-- RNA (RawSequence "GAUGGAACUUGACUACGUAAAUU")
rna :: Iso' RNA DNA
rna = rnaDnaIso

-- | Helper to make using 'rnaDnaIso' when converting from 'RNA' to 'DNA'.
--
-- >>> dna # RNA ("GAUGGAACUUGACUACGUAAAUU" ^. _RawSequence)
-- DNA (RawSequence "GATGGAACTTGACTACGTAAATT")
dna :: Iso' DNA RNA
dna = from rnaDnaIso
