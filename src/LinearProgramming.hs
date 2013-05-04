
module LinearProgramming
 (
 ) where

import LinearAlgebra
import qualified Data.Vector as V

m = fromLists [[1,0,1,-1],[0,1,-1,2]]

b = V.fromList [2,1] :: Vect

{- Solves the following
 - \[ \min c^T x \]
 - \[ \text{s.t.} Ax = b \]
 - \[ x >= 0 \]
 - There must be less or equal constraints than decision variables.
 - Dose not add slack variables
 - Returns Nothing if infeasible
 -}
simplexMethod :: Vect -> Matrix -> Vect -> Maybe Vect
simplexMethod c a b
    | V.length c /= (V.length . V.head) a || V.length a /= V.length b
                = error "Matrix and vector sizes do not match"
    | otherwise = Nothing

-- list of tuples of decision variable indices and their corresponding matrix
feasBases :: Matrix -> Vect -> [([Int], Matrix)]
feasBases mat b = filter (\(_,m) -> isFeasBasis (LinSys m b)) allBases
    where dvc       = (V.length $ V.head mat)
          conc      = V.length mat
          aInds     = combinations [0..(dvc-1)] conc
          allBases  = map (\ind -> (ind, colsFromList mat ind)) aInds

isFeasBasis :: LinearSystem -> Bool
isFeasBasis lns = V.and (V.map (>=0) (solLinSys lns))

-- Combination Related things
combinations xs 1   = map (\y -> [y]) xs
combinations (x:xs) i  = flatten $ map ((comb' i)) sl
    where sl                = subLists (x:xs)
          comb' j (y:ys)    = map (y :) (combinations ys (j-1))

flatten []      = []
flatten (x:xs)  = x ++ flatten xs

subLists []     = []
subLists [x]    = [[x]]
subLists (x:xs) = [(x:xs)] ++ subLists xs

