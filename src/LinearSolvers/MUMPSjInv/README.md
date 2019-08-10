# MUMPS

Efficiently solving large and sparse linear systems is essential in many areas of scientific computing. This project provides an interface to [MUMPS ( a MUltifrontal Massively Parallel sparse direct Solver)](http://mumps.enseeiht.fr/) in Julia. To improve the memory-efficiency of the method the software package  [METIS (Serial Graph Partitioning and Fill-reducing Matrix Ordering)](http://glaros.dtc.umn.edu/gkhome/views/metis) is used. Please refer to the MUMPS and metis websites for more information and updated versions.

## Copyright

Please make sure you have read and understood the licenses of [MUMPS  4.10.0](http://graal.ens-lyon.fr/MUMPS/index.php?page=dwnld)  and [Metis (4.0.3)](http://www.filewatcher.com/m/metis-4.0.3.tar.gz.522624-0.html). Due to the use of Metis, this code can be freely used for educational and research purposes by non-profit institutions and US government agencies only. The wrapper itself is provided as free software under MIT License in the hope it is useful.

Note that this software is provided "as is" with absolutely no warranty. Use at your own risk.


## Benchmark Example

In our experience MUMPS is considerably faster than Julia's backslash both in  solving real or complex systems. It also seems to consume less memory. Another benefit is that once the factorization is generated it can be applied to multiple right-hand sides. Here are results of two benchmarks ran on a MacBookPro with 2.6 GHz Intel Core i7 Processor and 16 GB RAM.

### Solving Poisson System
We are solving linear systems arising from discretization of Poisson's equation for different discretization sizes for a random right hand side. We compare the runtime for MUMPS solver and the backslash. 

| Grid size | MUMPS (sec) | Julia (sec) | Speedup|
| --------- | ------------|-------------|--------|
| 8^3  		|   0.0013    | 0.0030      | 2.243  |
| 16^3 		|   0.0197    | 0.0214      | 1.080  |
| 24^3 		|   0.1283    | 0.1112      | 0.866  |
| 32^3 		|   0.4160    | 0.3997      | 0.960  |
| 48^3 		|   2.3682    | 3.4196      | 1.443  |
| 64^3 		|   10.9231   | 16.1198     | 1.475  |

### Solving a complex system
As above, we obtain  discretizations of Poisson's equation for different grid sizes,  add an imaginary identity matrix and solve for one complex random right hand side. We compare the runtime for MUMPS solver and the backslash. 

| Grid size | MUMPS (sec) | Julia (sec) | Speedup|
|----------:|------------:|------------:|-------:|
| 8^3  		| 0.002       |   0.004     | 1.829  |
| 16^3 		| 0.047       |   0.087     | 1.867  |
| 24^3 		| 0.319       |   0.415     | 1.299  |
| 32^3 		| 1.312       |   2.502     | 1.907  |
| 48^3 		|11.091       |  24.057     | 2.168  |
| 64^3 		|58.472       | 157.315     | 2.690  |

Some more tests are provided in MUMPS/tests. 

## Installing MUMPS

This package was tested on Mac and Ubuntu only. Here, open a terminal navigate to MUMPS/src and run make. Please do a pull-request if you have improvements. 


**ToDo: Make MUMPS available from Julia's package manager! **





