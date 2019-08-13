[![Build Status](https://travis-ci.org/JuliaInv/jInv.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/jInv.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/jInv.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/jInv.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/0pxgtmm08b0w6wgh?svg=true)](https://ci.appveyor.com/project/JuliaInv/jinv-jl-81lel)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaInv.github.io/jInv.jl/latest)

# jInv

`jInv` is a flexible framework for PDE parameter estimation in Julia. It provides easy to extend core functions used in PDE-constrained inverse problems.
Our goal is to solve parameter estimation problems efficiently and in parallel. For more details see (http://arxiv.org/abs/1606.07399).


# Overview

jInv consists of five submodules:

1. `ForwardShare` - methods for solving forward problems in parallel.
2. `InverseSolve` - methods commonly used in inverse problems such as misfit functions, regularization and numerical optimization.
3. `Mesh` - regular and tensor meshes in 2D and 3D as well as differential operators.
4. `LinearSolvers` - interfaces to sparse and (if installed) direct linear solvers that can be used for solving the discretized PDEs.
5. `Utils` - utility functions

# Requirements

jInv is intended for use with Julia versions 0.7. and requires:

1. [`KrylovMethods.jl`](https://github.com/lruthotto/KrylovMethods.jl)  - iterative methods for solving (sparse) linear systems.

Additional (optional) packages for higher performance. 
[`ParSpMatVec.jl`](https://github.com/lruthotto/ParSpMatVec.jl) - shared memory implementation for sparse matrix vector products. `jInv` detects automatically if this package is properly set and uses it by default.

Additional add-on package:

1. [`MUMPSjInv.jl`](https://github.com/JuliaInv/MUMPSjInv.jl) - wrapper for MUMPS. Used as a direct linear solver.
2. ['Pardiso.jl'](https://github.com/JuliaSparse/Pardiso.jl) - The extension using this package exists but was not tested for Julia version 1.0.

The mesh module in jInv features regular and tensor meshes but can also be extended by
1. [`JOcTree`](https://github.com/JuliaInv/JOcTree) - efficient spatially adaptive meshes

# Installation

In julia type:
```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.test("jInv")
```

# Examples

Some inversion examples can be found in the `examples` folder. 

# Packages using jInv

jInv provides generic methods for PDE parameter estimation problems. In order to use it for applications, some methods need to be extended by specifying solvers of the forward problem, sensitivity matrix vector products, etc. This should be organized in small modules. Currently implemented are the following PDE models:

1. [`DivSigGrad.jl`](https://github.com/JuliaInv/DivSigGrad.jl) - Inverse conductivity problems in statics
2. [`jInvSeismic.jl`](https://github.com/JuliaInv/jInvSeismic.jl) - Seismic inversion packages: Full Waveform Inversion (jInvSeismic.FWI) and travel time tomography (jInvSeismic.EikonalInv)
3. [`MaxwellFrequency`](https://github.com/JuliaInv/MaxwellFrequency) - Inversion for conductivity in Maxwell's equations

# Acknowledgements

This material is in part based upon work supported by the National Science Foundation under Grant Number 1522599. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
