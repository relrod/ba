{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types.RawSequence
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- Types which represent a raw sequence of some kind.
----------------------------------------------------------------------------
module Bio.Algorithm.Types.RawSequence where

import Control.Lens
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Text.Lazy as TL
import qualified Data.Text.Lazy.Encoding as TL
import Data.String

newtype RawSequence = RawSequence BL.ByteString deriving (Eq, Ord, Show)

instance IsString RawSequence where
  fromString = RawSequence . BL.pack
  {-# INLINE fromString #-}

instance Reversing RawSequence where
  reversing (RawSequence s) = RawSequence . reversing $ s
  {-# INLINE reversing #-}

class AsRawSequence s where
  _RawSequence :: Iso' s RawSequence

instance AsRawSequence RawSequence where
  _RawSequence = id
  {-# INLINE _RawSequence #-}

instance AsRawSequence BL.ByteString where
  _RawSequence = iso RawSequence (\(RawSequence r) -> r)
  {-# INLINE _RawSequence #-}

instance AsRawSequence String where
  _RawSequence = iso (RawSequence . BL.pack) (\(RawSequence r) -> BL.unpack r)
  {-# INLINE _RawSequence #-}

-- | NOTE: This uses 'TL.decodeUtf8' and 'TL.encodeUtf8'.
instance AsRawSequence TL.Text where
  _RawSequence = iso (RawSequence . TL.encodeUtf8) (\(RawSequence r) -> TL.decodeUtf8 r)
  {-# INLINE _RawSequence #-}
