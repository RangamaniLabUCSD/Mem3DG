# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)

[![Build Status](https://travis-ci.com/RangamaniLabUCSD/DDG_membrane.svg?token=HxusyqZoDyxkhvY6GCzF&branch=master)](https://travis-ci.com/RangamaniLabUCSD/DDG_membrane)

3-D computational model for lipid membrane 

## Acknowledging your use of Mem3DG
Thanks for using Mem3DG! The developers would love to hear how you are using the tool. Please send us an email or post on GitHub letting us know.

Please cite:

## Installation

```
mkdir build
cd build
cmake -DBUILD_PYDDG=ON -DWITH_NETCDF=ON -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

## Dependencies

* Optional trajectory output uses NetCDF. Both the NetCDF and NetCDF_cxx4 libraries must be available on your system to use this feature.
