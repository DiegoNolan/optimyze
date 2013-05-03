
import LinearAlgebra
import System.Random
import qualified Data.Vector as V


testSolveLinSys dim = do
    putStrLn "Solving the following Linear System"
    (LinSys m b) <- randomLinSys dim
    (putStr . show) (LinSys m b)
    putStrLn "Solution"
    let sol = (solLinSys (LinSys m b))
    (putStrLn . show) sol
    putStrLn "Matrix times solution"
    (putStrLn . show) (multMV m sol)

randomLinSys :: Int -> IO LinearSystem
randomLinSys s = do
    mat <- randomMatrix s s
    vec <- randomVect s
    return (LinSys mat vec)

randomMatrix :: Int -> Int -> IO Matrix
randomMatrix _ 0 = return V.empty
randomMatrix r c = do
    row <- randomVect r
    rest <- randomMatrix r (c-1)
    let mat = row `V.cons` rest
    return mat

randomVect :: Int -> IO Vect
randomVect dim = do
    gen <- newStdGen
    return $ V.fromList (map fromIntegral  $ take dim $
                randomRs (-10::Int,10::Int) gen)

