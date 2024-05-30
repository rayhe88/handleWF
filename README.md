# handleWF
A code to manipulate the wave function of an electronic structure calculation in WFX format.


This is the version 2 written using C++.

## Compilation

Create the `build` subdirectory
```
mkdir build
```

Move to the `build`
```
cd build
```

Type the next commands
```
cmake -H. -DUSE_SYCL=On -DCMAKE_CXX_COMPILER=icpx ../src/
```

Finally build the project
```
cmake --build .
```

## Testing

```
./handleWF.x ../test/dimer_HCOOH.wfx
Version: 0.0
Compilation Date: Apr  9 2024  13:37:02
Git SHA1: 6dba7e6
 Points ( 80,80,40)
 TotalPoints : 256000
 Running on 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
 Points ( 80,80,40)
 TotalPoints : 256000
 Running on 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
 Points ( 80,80,40)
 TotalPoints : 256000
 Time for CPU : 5.0914e+07 μs
 Time for GPU  : 1.39384e+06 μs (Kernel 1)
 Time for GPU  : 2.16659e+06 μs (Kernel 2)
 Ratio CPU/GPU (kernel1) : 36.5278
 Ratio CPU/GPU (kernel2) : 23.4997

```

## Acknowledgements
This research used resources of the Argonne Leadership Computing Facility, which is a DOE Office of Science User Facility supported under Contract DE-AC02-06CH11357. Argonne National Laboratory’s work was supported by the U.S. Department of Energy, Office of Science, under contract DE-AC02-06CH11357.


All rights reserved. Copyright Argonne National Laboratory UChicago LLC. Raymundo Hernandez-Esparza

