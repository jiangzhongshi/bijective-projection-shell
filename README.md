# Bijective Projection in a Shell


Zhongshi Jiang, Teseo Schneider, Denis Zorin, Daniele Panozzo. 
*ACM Transactions on Graphics (In Proceedings of SIGGRAPH Asia 2020)*

<img src="https://i.imgur.com/sgiVMlh.jpg" width="200"/>

## Overview
We introduce an algorithm to convert a self-intersection free, orientable, and manifold triangle mesh T into a generalized prismatic shell equipped with a bijective projection operator to map T to a class of discrete surfaces contained within the shell whose normals satisfy a simple local condition. Properties can be robustly and efficiently transferred between these surfaces using the prismatic layer as a common parametrization domain.

The combination of the prismatic shell construction and corresponding projection operator is a robust building block readily usable in many downstream applications, including the solution of PDEs, displacement maps synthesis, Boolean operations, tetrahedral meshing, geometric textures, and nested cages.

## Installation
Prerequisites: `cmake, boost, hdf5, gmp, mpfr`
Also refer to [cmake.yml](.github/workflows/cmake.yml).
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

## Usage
We use [HDF5](https://www.hdfgroup.org/solutions/hdf5/) for serialize results. And the fields `mtop` (top surface), `mV` (midsurface), `mbase` (bottom surface), `mF` (connectivity) represent the output shell.

Some scripts are provided in [FigureScripts](FigureScripts.md).
Please stay tuned for more information!

