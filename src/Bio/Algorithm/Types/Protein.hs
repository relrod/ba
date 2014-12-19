{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}

-----------------------------------------------------------------------------
-- |
-- Module : Bio.Algorithm.Types.Protein
-- Copyright : (C) 2014 Ricky Elrod
-- License : BSD2 (see LICENSE file)
-- Maintainer : Ricky Elrod <ricky@elrod.me>
-- Stability : experimental
-- Portability : lens
--
-- Protein types
----------------------------------------------------------------------------
module Bio.Algorithm.Types.Protein where

--import Control.Lens
--import qualified Data.ByteString.Lazy.Char8 as BL
--import qualified Data.Text.Lazy as TL
--import qualified Data.Text.Lazy.Encoding as TL

data Protein = Protein Char | Stop deriving (Eq, Ord, Show)
