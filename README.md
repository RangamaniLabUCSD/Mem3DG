# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)

[![Build Status](https://travis-ci.com/RangamaniLabUCSD/Mem3DG.svg?token=HxusyqZoDyxkhvY6GCzF&branch=master)](https://travis-ci.com/RangamaniLabUCSD/Mem3DG)

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

## Temporary notes for setting up netcdf (especially on windows...)

1. Download `vcpkg` and follow the instructions to install
2. Install 32 or 64 bit version of `netcdf-c` and `netcdf-cxx4` libraries depending upon your configuration.

   `vcpkg install netcdf-c:x64-windows netcdf-cxx4:x64-windows`

   Remove the `:x64-windows` from the above string for the 32 bit libraries.

3. Configure the vcpkg CMake toolchain `vcpkg integrate install`
4. Copy and paste the `-DCMAKE_TOOLCHAIN...` string as an input into your CMake configuration.
5. Build as normal

1.  python setup.py build -- -DCMAKE_TOOLCHAIN_FILE="C:/Users/Kieran/vcpkg/scripts/buildsystems/vcpkg.cmake" -G "Visual Studio 16 2019" -T host=x64 -A x64 -- /m:6

## Dependencies

* Optional trajectory output uses NetCDF. Both the NetCDF and NetCDF_cxx4 libraries must be available on your system to use this feature.
