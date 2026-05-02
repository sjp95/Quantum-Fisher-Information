#g++ -std=c++11 st1.cpp -I /usr/local/opt/lapack/include -L/usr/local/opt/lapack/lib -llapack -o h -fopenmp
#g++ -std=c++11 st1.cpp -llapack -o h -fopenmp

g++ -std=c++11 st1.cpp -o h -llapack -I /usr/include/lapacke \-L /usr/lib64  -fopenmp