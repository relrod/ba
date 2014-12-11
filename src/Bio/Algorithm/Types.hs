{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Bio.Algorithm.Types where

import Control.Applicative
import Control.Lens
import qualified Data.ByteString.Char8 as B
import qualified Data.Text as T
import qualified Data.Text.Encoding as T

newtype RawSequence = RawSequence B.ByteString deriving (Eq, Ord, Show)

class AsRawSequence p f s where
  _RawSequence :: Optic' p f s RawSequence

instance AsRawSequence p f RawSequence where
  _RawSequence = id

instance (Choice p, Applicative f) => AsRawSequence p f B.ByteString where
  _RawSequence = iso RawSequence (\(RawSequence r) -> r)

instance (Choice p, Applicative f) => AsRawSequence p f String where
  _RawSequence = iso (RawSequence . B.pack) (\(RawSequence r) -> B.unpack r)

-- | NOTE: This uses 'T.decodeUtf8' and 'T.encodeUtf8'.
instance (Choice p, Applicative f) => AsRawSequence p f T.Text where
  _RawSequence = iso (RawSequence . T.encodeUtf8) (\(RawSequence r) -> T.decodeUtf8 r)
