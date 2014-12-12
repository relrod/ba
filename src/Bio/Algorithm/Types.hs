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
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Text.Lazy as TL
import qualified Data.Text.Lazy.Encoding as TL

newtype RawSequence = RawSequence BL.ByteString deriving (Eq, Ord, Show)

class AsRawSequence p f s where
  _RawSequence :: Optic' p f s RawSequence

instance AsRawSequence p f RawSequence where
  _RawSequence = id
  {-# INLINE _RawSequence #-}

instance (Profunctor p, Functor f) => AsRawSequence p f BL.ByteString where
  _RawSequence = iso RawSequence (\(RawSequence r) -> r)
  {-# INLINE _RawSequence #-}

instance (Profunctor p, Functor f) => AsRawSequence p f String where
  _RawSequence = iso (RawSequence . BL.pack) (\(RawSequence r) -> BL.unpack r)
  {-# INLINE _RawSequence #-}

-- | NOTE: This uses 'TL.decodeUtf8' and 'TL.encodeUtf8'.
instance (Profunctor p, Functor f) => AsRawSequence p f TL.Text where
  _RawSequence = iso (RawSequence . TL.encodeUtf8) (\(RawSequence r) -> TL.decodeUtf8 r)
  {-# INLINE _RawSequence #-}
