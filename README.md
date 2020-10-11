# Introduction:
The first task for the Parallel Computing course in CMC MSU. Implement the
conjugate gradient method. Parallelize it using OpenMP.
# Building a program:
Run a "tsk1" target on make.

To make the separate measurements for generate, fill and solver phases, run a
"tsk1\_Measure" target. It generates tsk1\_msr executable.

To make the separate measurement for the basic operations of the solver, run a
"tsk1\_Measure\_Solver" target. It generates tsk1\_msr\_slv executable.

# Launching a program
The first argument is an input file. The threads number is specified with "-t"
option. Enable the debug print with "-d" option.

# Perfomance results
I measured the perfomance on the Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz. It
has 4 cores.

A perfomance for the different phases.

|Measure |Matrix size |Threads|Time (s.) |
|--------|------------|-------|----------|
|All     |       50000|      1| 0.0800002|
|All     |       50000|      4| 0.0450001|
|All     |      500000|      1| 0.0800002|
|All     |      500000|      4| 0.0450001|
|All     |     5000000|      1| 7.027    |
|All     |     5000000|      4| 3.142    |
|Generate|       50000|      1| 0.0040019|
|Generate|       50000|      4| 0.0020009|
|Generate|      500000|      1| 0.0179999|
|Generate|      500000|      4| 0.0089998| 
|Generate|     5000000|      1| 0.201    |
|Generate|     5000000|      4| 0.0649998|
|Fill    |       50000|      1| 0.0150001|
|Fill    |       50000|      4| 0.0069997|
|Fill    |      500000|      1| 0.1530000|
|Fill    |      500000|      4| 0.0439999|
|Fill    |     5000000|      1| 1.558    |
|Fill    |     5000000|      4| 0.425    |
|Solver  |       50000|      1| 0.0540001|
|Solver  |       50000|      4| 0.0330000|
|Solver  |      500000|      1| 0.4590000|
|Solver  |      500000|      4| 0.2640000| 
|Solver  |     5000000|      1| 5.265    |
|Solver  |     5000000|      4| 2.65     | 


The memory usage for the different phases:

|Measure |Matrix size |Memory     |
|--------|------------|-----------|
|All     |       50000|    9408512|
|All     |      500000|   80035840| 
|All     |     5000000|  783544320| 
|Generate|       50000|    7503872| 
|Generate|      500000|   74911744| 
|Generate|     5000000|  742342656| 
|Fill    |       50000|     405504|
|Fill    |      500000|    4050944|
|Fill    |     5000000|   40128512| 
|Solver  |       50000|     417792|
|Solver  |      500000|          0| 
|Solver  |     5000000|          0| 

The operations inside the solve also can be measured separately
DotPr — the dot product
LinComb - the linear combination
SparseMV - the sparse multiplication

|Measure |Matrix size |Threads|Time (s.) |
|--------|------------|-------|----------|
|Solver  |       50000|      1| 0.0540001|
|Solver  |       50000|      4| 0.0330000|
|Solver  |      500000|      1| 0.4590000|
|Solver  |      500000|      4| 0.2640000| 
|Solver  |     5000000|      1| 5.265    |
|Solver  |     5000000|      4| 2.65     | 
|DotPr   |       50000|      1| 0.0039997|
|DotPr   |       50000|      4| 0.0019999|
|DotPr   |      500000|      1| 0.0280015|
|DotPr   |      500000|      4| 0.019001 |
|DotPr   |     5000000|      1| 0.345    |
|DotPr   |     5000000|      4| 0.161    |
|LinComb |       50000|      1| 0.0109994|
|LinComb |       50000|      4| 0.0069995|
|LinComb |      500000|      1| 0.0749989|
|LinComb |      500000|      4| 0.0400007|
|LinComb |     5000000|      1| 0.881    |
|LinComb |     5000000|      4| 0.443    |
|SparseMV|       50000|      1| 0.0250008|
|SparseMV|       50000|      4| 0.0170002|
|SparseMV|      500000|      1| 0.258    |
|SparseMV|      500000|      4| 0.0900004|
|SparseMV|     5000000|      1| 2.903    | 
|SparseMV|     5000000|      4| 1.072    |

The memory usage:

|Measure |Matrix size |Memory     |
|--------|------------|-----------|
|All     |       50000|    9408512|
|All     |      500000|   80035840| 
|Generate|       50000|    7503872| 
|Generate|      500000|   74911744| 
|Fill    |       50000|     405504|
|Fill    |      500000|    4050944|
|Solver  |       50000|     417792|
|Solver  |      500000|          0| 

TODO: Measure the perfomance on the cluster.
