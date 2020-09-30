#include <iostream>
#include <netcdf>

#include "mem3dg/solver/icosphere.h"
#include "mem3dg/solver/trajfile.h"
#include "mem3dg/solver/util.h"

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

// We are writing 2D data, a 6 x 12 grid
constexpr int nx = 6;
constexpr int ny = 12;

// Return this in event of a problem
constexpr int nc_err = 2;

int main() {
  namespace gc = ::geometrycentral;
  namespace gcs = ::geometrycentral::surface;

  std::vector<gc::Vector3> coords;
  std::vector<std::vector<std::size_t>> polygons;

  ddgsolver::tetrahedron(coords, polygons);

  gcs::SimplePolygonMesh soup(polygons, coords);
  soup.mergeIdenticalVertices();

  std::unique_ptr<gcs::SurfaceMesh> ptrMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> ptrVpg;
  std::tie(ptrMesh, ptrVpg) =
      gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);

  auto file = ddgsolver::TrajFile::newFile("test.nc", *ptrMesh, true);

  file.writeTime(file.getNextFrameIndex(), 1);
  file.writeTime(file.getNextFrameIndex(), 2);

  file.writeCoords(
      0, gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions));
  file.writeCoords(
      3, gc::EigenMap<double, 3>(ptrVpg->inputVertexPositions));

  double x, y;
  ddgsolver::TrajFile::EigenVector vec1, vec2;

  std::tie(x, vec1) = file.getTimeAndCoords(0);
  std::cout << "Time " << x << std::endl << vec1 << std::endl;
  
  auto file2 = ddgsolver::TrajFile::openReadOnly("test.nc");
  std::tie(y, vec2) = file2.getTimeAndCoords(1);
  std::cout << "Time " << y << std::endl << vec2 << std::endl;

  std::cout << "EOF" << std::endl;
  return 0;

  // The default behavior of the C++ API is to throw an exception if
  // an error occurs
  try {
    // This is the data array we will write. It will just be filled
    // with a progression of numbers for this example.
    int dataOut[nx][ny];

    // Create some pretend data. If this wasn't an example program, we
    // would have some real data to write, for example, model output.
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        dataOut[i][j] = i * ny + j;
      }
    }

    // Create the file. The Replace parameter tells netCDF to overwrite
    // this file, if it already exists.
    netCDF::NcFile dataFile("simple_xy.nc", netCDF::NcFile::replace);

    // Create netCDF dimensions
    auto framesDim = dataFile.addDim("frame");
    auto xDim = dataFile.addDim("x", nx);
    auto yDim = dataFile.addDim("y", ny);

    // Define the variable. The type of the variable in this case is
    // ncInt (32-bit integer)
    auto data =
        dataFile.addVar("coordinates", netCDF::ncInt, {framesDim, xDim, yDim});

    dataFile.putAtt("TestAttribute", "FOO");

    data.putAtt("units", "angstroms");

    // Write the data to the file. Although netCDF supports reading
    // and writing subsets of data, in this case we write all the data
    // in one operation.
    data.putVar({0, 0, 0}, {1, 6, 12}, &dataOut);
    data.putVar({2, 0, 0}, {1, 6, 12}, &dataOut);
    
    // The file will be automatically close when the NcFile object goes
    // out of scope. This frees up any internal netCDF resources
    // associated with the file, and flushes any buffers.
  } catch (netCDF::exceptions::NcException &e) {
    std::cout << e.what() << std::endl;
    return nc_err;
  }

  // Now read the data back in
  try {
    // This is the array we will read into
    int dataIn[nx][ny];

    // Open the file for read access
    netCDF::NcFile dataFile("simple_xy.nc", netCDF::NcFile::read);

    // Retrieve the variable named "data"
    auto data = dataFile.getVar("coordinates");
    if (data.isNull()){
      std::cout << "Null data" << std::endl;
      return nc_err;
    }
    data.getVar(dataIn);

    // Check the values.
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        if (dataIn[i][j] != i * ny + j) {
          std::cout << "Data mismatch" << std::endl;
          return nc_err;
        }
      }
    }
  } catch (netCDF::exceptions::NcException &e) {
    std::cout << e.what() << std::endl;
    return nc_err;
  }
}
