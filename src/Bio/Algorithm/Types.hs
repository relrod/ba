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
-- >>> T.pack "CCTTGGAA" ^. lazy . _RawSequence . to (dnaReverseComplement . DNA)
-- DNA (RawSequence "TTCCAAGG")
--
-- Note that you can convert between lazy and strict versions of ByteString and
-- Text by using the appropriate
-- <https://hackage.haskell.org/package/lens/docs/Control-Lens-Iso.html#t:Strict lens combinators>.
--
-- This module simply re-exports all of Bio.Algorithm.Types.*.
----------------------------------------------------------------------------
module Bio.Algorithm.Types (
  module Bio.Algorithm.Types.DNA
, module Bio.Algorithm.Types.RawSequence
, module Bio.Algorithm.Types.RNA
) where

import Bio.Algorithm.Types.DNA
import Bio.Algorithm.Types.RawSequence
import Bio.Algorithm.Types.RNA

-- $setup
-- >>> import Bio.Algorithm.Sequence
-- >>> import Control.Lens
-- >>> import qualified Data.Text as T
