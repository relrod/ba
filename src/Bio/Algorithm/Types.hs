{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- An experimental approach to working with various bioinformatics types in
-- Haskell. Makes very heavy use of lens.
--
-- Making such strong use of lenses, we enable an API that lets you do things
-- like this:
--
-- >>> T.pack "CCTTGGAA" ^. _RawSequence . to reverseComplement
-- RawSequence "TTCCAAGG"
-- it :: RawSequence
----------------------------------------------------------------------------
module Bio.Algorithm.Types where

import Control.Lens
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Text.Encoding as T

newtype RawSequence = RawSequence B.ByteString deriving (Eq, Ord, Show)

class AsRawSequence p f s where
  _RawSequence :: Optic' p f s RawSequence

instance AsRawSequence p f RawSequence where
  _RawSequence = id

instance (Profunctor p, Functor f) => AsRawSequence p f B.ByteString where
  _RawSequence = iso RawSequence (\(RawSequence r) -> r)

instance (Profunctor p, Functor f) => AsRawSequence p f String where
  _RawSequence = iso (RawSequence . B.pack) (\(RawSequence r) -> B.unpack r)

-- | NOTE: This uses 'T.decodeUtf8' and 'T.encodeUtf8'.
instance (Profunctor p, Functor f) => AsRawSequence p f T.Text where
  _RawSequence = iso (RawSequence . T.encodeUtf8) (\(RawSequence r) -> T.decodeUtf8 r)
