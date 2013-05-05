
module LinearProgramming
 (
 ) where

import LinearAlgebra
import qualified Data.Vector as V
import Data.List (delete)

ctest = V.fromList [-1,0,0,2] :: Vect

mtest = fromLists [[1,0,1,-1],[0,1,-1,2]]

btest = V.fromList [2,1] :: Vect

ltest = LinSys mtest btest

{- Solves the following
 - \[ \min c^T x \]
 - \[ \text{s.t.} Ax = b \]
 - \[ x >= 0 \]
 - There must be less or equal constraints than decision variables.
 - Dose not add slack variables
 - Returns Nothing if infeasible
 -}
--simplexMethod :: Vect -> Matrix -> Vect -> Maybe Vect
simplexMethod c a b
    | V.length c /= (V.length . V.head) a || V.length a /= V.length b
                = error "Matrix and vector sizes do not match"
    | otherwise = simplexMeth c a b (allBases a b)

--simplexMeth :: Vect -> Matrix -> Vect -> [([Int],Matrix)] -> Maybe Vect
simplexMeth _ _ _ []    = Nothing -- No feasible basis
simplexMeth c a b bases
    | (not . V.and) $ V.map (>=0) xb = simplexMeth c a b (drop 1 bases)
    | V.and $ V.map (<=0) rc         = Just (fullVector xb xbi xni)
    | otherwise                      = simplexMeth c a b ((newbasi,newB): drop 1 bases)
        where (xbi,bas) = head bases
              xb        = solLinSys (LinSys bas b)
              xni       = indexComp xbi ((V.length $ V.head a) - 1)
              nonb      = transpose $ colsFromList a xni
              cb        = subFromList c xbi
              cn        = subFromList c xni 
              mults     = solLinSys $ LinSys (transpose bas) cb -- Lagrange mults
              rc        = cn `minusVV` (multMV nonb mults)
              ent       = indOfMax rc -- index to enter basis
              d         = solLinSys $ LinSys bas (column a ent) -- direction
              ratio     = V.zipWith (/) xb d
              leave     = xbi !! (indOfMax ratio)
              newbasi   = addOrd (delete leave xbi) ent
              newB      = colsFromList a newbasi

fullVector :: Vect -> [Int] -> [Int] -> Vect
fullVector vs [] ns     = V.replicate (length ns) 0
fullVector vs bs []     = vs
fullVector vs (b:bs) (n:ns)
    | b < n             = V.head vs `V.cons` fullVector (V.drop 1 vs) (bs) (n:ns)
    | otherwise         = 0 `V.cons` fullVector vs (b:bs) ns

addOrd :: [Int] -> Int -> [Int]
addOrd [] a      = [a]
addOrd [x]  a    = if x > a then a : [x] else x : [a]
addOrd (x:y:ys) a
    | a > x && a < y    = x : a : y : ys
    | otherwise         = addOrd (y:ys) a

-- takes a list of indexs and an upper bound and returns the lists compliment
-- [1,2,4] 5 = [0,3,5]
indexComp :: [Int] -> Int -> [Int]
indexComp ind up = indexComp' ind up 0
    where indexComp' []  up j   =  [j..up]
          indexComp' (i:is) up j
                | j > up        = []
                | j == i        = indexComp' is up (j+1)
                | otherwise     = j : indexComp' (i:is) up (j+1)

-- list of tuples of decision variable indices and their corresponding matrix
allBases :: Matrix -> Vect -> [([Int], Matrix)]
allBases mat b = map (\ind -> (ind, colsFromList mat ind)) aInds
    where dvc       = (V.length $ V.head mat)
          conc      = V.length mat
          aInds     = combinations [0..(dvc-1)] conc

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

