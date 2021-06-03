---
layout: splash
title: Mem3DG 
hide-name: true
---

[![Build Status](https://travis-ci.com/RangamaniLabUCSD/Mem3DG.svg?token=HxusyqZoDyxkhvY6GCzF&branch=master)](https://travis-ci.com/RangamaniLabUCSD/Mem3DG)
[![PyPI](https://img.shields.io/pypi/v/pymem3dg)](https://pypi.org/project/pymem3dg/)

Mem3DG is a flexible software package to model the membrane and its dynamics using unstructured meshes.
This work is currently under heavy development, please star this repository to follow along!

## Acknowledging your use of Mem3DG

Mem3DG is developed by [Cuncheng Zhu](https://github.com/cuzhucuncheng), [Christopher T. Lee](https://ctlee.github.io/), with contributions from others.
Development of Mem3DG is funded in part by AFOSR MURI FA9550-18-1-0051, and a Hartwell Foundation Postdoctoral Fellowship.

## Installation

```
git submodule --init --recursive
mkdir build
cd build
cmake -DBUILD_PYDDG=ON -DWITH_NETCDF=ON -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

Source released can also be obtained from [PyPi](https://pypi.org/project/pymem3dg/).

## Dependencies

* Meshes are represented using [Geometry-Central](https://geometry-central.net/).
* Optional trajectory output uses [NetCDF-cxx4](https://github.com/Unidata/netcdf-cxx4).
