{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types.DNA
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- DNA sequences
----------------------------------------------------------------------------
module Bio.Algorithm.Types.DNA where

import Bio.Algorithm.Types.RawSequence
import Control.Lens

newtype DNA = DNA RawSequence deriving (Eq, Ord, Show)

instance (Profunctor p, Functor f) => AsRawSequence p f DNA where
  _RawSequence = iso (\(DNA r) -> r) DNA
  {-# INLINE _RawSequence #-}
