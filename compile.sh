find . -name "*~" -type f -delete
make test-debug -j 4
./bin/debug/DistMatrix 1 1 16 16 
./bin/debug/BLAS/Gemm 1 1 N N 128 128 128 32 1 1 1 0
./bin/debug/BLAS/Gemv 1 1 8 2 1
./bin/debug/BLAS/Ger 1 1 8 2 1
