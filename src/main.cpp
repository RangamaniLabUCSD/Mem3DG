
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <iostream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "ddgsolver/ddg_solver.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/force.h"


namespace gc  = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// overload << to print vector;
template<typename T>
std::ostream& operator<< (std::ostream& output, const std::vector<T>& v) {
	output << "[";
	for (size_t i = 0; i != v.size() - 1; ++i) {
		output << v[i] << ",";
	}
	output << v[v.size()-1];
	output  << "]";
	return output;

}

int main() {

	// Construct an icosphere
	std::vector<gc::Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons;
	icosphere(coords, polygons, 1);

	// Construct a mesh
	gc::PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();

	// halfedge mesh and vertex postion geometry from face and vertex coords 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

	// initiate force object f
	force f(*ptrmesh, *ptrvpg);
	std::cout << "Sizeof Force: " << sizeof(f) << std::endl;

	// calculate the bending force 
	Eigen::Matrix<double, Eigen::Dynamic, 3> bf = f.bending_force(1.0, 0);
	//gcs::IntrinsicGeometryInterface& geometry = *ptrvpg;

	// calculate the stretching force
	Eigen::Matrix<double, Eigen::Dynamic, 3> sf = f.stretching_force(1.0, 1.0);

	 /*polyscope::init();
	 polyscope::registerSurfacemesh("myMesh",
	 							   ptrvpg->inputVertexPositions,
	 							   ptrmesh->getFaceVertexList());
	 polyscope::show();*/

	return 0;
}


