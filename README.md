# TSP using MPI
This project implements MPI for the travelling salesperson problem. It explores distributed scheduling and load balancing the calculation as well as expirementing with efficient calculations.

To build
```Bash
mkdir build
cd build
cmake ..
make -j 
```

To run
```Bash
cd build
mpirun -n ${number of cores/computers} TSP_MPI 14 ../cities/city14.txt
```
