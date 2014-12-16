{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types.RNA
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- RNA sequences
----------------------------------------------------------------------------
module Bio.Algorithm.Types.RNA where

import Bio.Algorithm.Types.RawSequence
import Control.Lens

newtype RNA = RNA RawSequence deriving (Eq, Ord, Show)

instance AsRawSequence RNA where
  _RawSequence = iso (\(RNA r) -> r) RNA
  {-# INLINE _RawSequence #-}
