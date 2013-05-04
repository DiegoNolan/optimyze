
module LinearAlgebra
 ( Vect
 , Matrix
 , LinearSystem (LinSys)
 , fromLists
 , innerProduct
 , transpose
 , multVM
 , multMV
 , multMM
 , scalarM
 , scalarV
 , column
 , colsFromList
 , product
 , solLinSys
 ) where

import qualified Data.Vector as V

-- Assumes the linear system has only one solution, ie
-- square matrix and linear independence, matrix and vector
-- size must match

type Vect = V.Vector Double

type Matrix = V.Vector Vect

data LinearSystem = LinSys Matrix Vect

instance Show LinearSystem where
    show (LinSys mat v) = allRows mat v
        where allRows mat v = foldl (\str x -> str ++ rowToString mat v x) []
                    [0..((V.length v)-1)]
              rowToString mat v r = (show $ V.toList (mat V.! r)) ++ " " ++
                    (show (v V.! r)) ++ "\n"

fromLists :: [[Double]] -> Matrix
fromLists vs = V.fromList $ map V.fromList vs

-- General
innerProduct :: Vect -> Vect -> Double
innerProduct v1 v2 | V.length v1 /= V.length v2 =
    error "Not same length vectors"
innerProduct v1 v2 = V.sum $ V.zipWith (*) v1 v2

transpose :: Matrix -> Matrix
transpose m = V.map (column m) (V.fromList [0..(V.length (m V.! 0) - 1)])

multVM :: Vect -> Matrix -> Vect
multVM v m = V.map (innerProduct v) (transpose m)

multMV :: Matrix -> Vect -> Vect
multMV m v = V.map (innerProduct v) m

multMM :: Matrix -> Matrix -> Matrix
multMM f s = V.map (\r -> multVM r s) f

scalarM :: Double -> Matrix -> Matrix
scalarM s m = V.map (scalarV s) m

scalarV :: Double -> Vect -> Vect
scalarV s v = V.map (*s) v

column :: Matrix -> Int -> Vect
column mat j | j < 0 || j >= V.length (mat V.! 0) = error "Out of range"
column mat j = column' mat j 0
    where column' mat j i
                | i >= V.length mat = V.empty
                | otherwise         = ((mat V.! i) V.!j)
                        `V.cons` (column' mat j (i+1))

colsFromList :: Matrix -> [Int] -> Matrix
colsFromList m inds = transpose ( V.fromList $ map (column m) inds)

-- Solving Linear Systems of Equations
solLinSys :: LinearSystem -> Vect
solLinSys (LinSys mat b) | V.length mat /= V.length (mat V.! 0) = error
    "Has to be a square matrix"
solLinSys (LinSys mat b) | V.length (mat V.! 0) /= V.length b = error
    "Matrix and Vector sizes don't match"
solLinSys lns = (solEch . toEchelon) lns
    where solEch lns = sol' lns V.empty

-- Solves Linear System in Echelon form
sol' (LinSys m b) s | V.length s == V.length b  = s
sol' (LinSys m b) s = sol' (LinSys m b) (c `V.cons` s)
    where i = (V.length b) - (V.length s) - 1
          r = V.drop (i+1) (m V.! i)
          c = ((b V.! i) - (innerProduct s r))/((m V.! i) V.! i)

-- Converts a linear system to Echelon orm
toEchelon :: LinearSystem -> LinearSystem
toEchelon (LinSys mat b)
    | V.length b == 1   = LinSys mat b
    | otherwise         = LinSys mm bb
        where mm = (V.head mat) `V.cons` (V.map (0 `V.cons`) mret)
              bb = (V.head b) `V.cons` bret
              (LinSys mret bret) = ((toEchelon . lowerRight)
                    ( elimCols (LinSys mat b) 1))

lowerRight :: LinearSystem -> LinearSystem
lowerRight (LinSys m b) = LinSys (V.map (V.drop 1) ((V.drop 1) m)) (V.drop 1 b)

-- turns the front colum of the matrix in a linear system to zero
-- returns the equivalent linear system
elimCols :: LinearSystem -> Int -> LinearSystem
elimCols (LinSys m b) r
    | r >= V.length m       = LinSys m b
    | otherwise             = elimCols (elimCol (LinSys m b) r) (r+1)

elimCol :: LinearSystem -> Int -> LinearSystem
elimCol (LinSys mat b) r = LinSys mret vret
    where m         = (-(V.head $ mat V.! r)/(V.head $ mat V.! 0))
          mbegin    = V.take r mat
          mrest     = V.drop (r+1) mat
          mret      = mbegin V.++ row `V.cons` mrest
          vbegin    = V.take r b
          vrest     = V.drop (r+1) b
          vret      = vbegin V.++ ((b V.! 0)*m + (b V.! r)) `V.cons` vrest
          cnclRow   = V.map (*m) (mat V.! 0)
          row       = V.zipWith (+) cnclRow (mat V.! r)

