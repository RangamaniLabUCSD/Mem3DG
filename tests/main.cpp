
#include <iostream>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "ddgsolver/icosphere.h"
#include "ddgsolver/force.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"


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

	//initialization of the geometry and mesh
	// Construct an icosahedron
	// Construct a std::vector of Vector 3
	std::vector<gc::Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons;

	icosphere(coords, polygons, 3);

	gc::PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();

	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

	auto& mesh = *ptrmesh;
	auto& vpg  = *ptrvpg;

    auto pos = mapVecToEigen<double, 3>(vpg.inputVertexPositions.rawdata());

    std::cout << "decltype(pos) is " << type_name<decltype(pos)>() << '\n';

    for (std::size_t i = 0; i < mesh.nVertices(); ++i){
    	std::cout << vpg.inputVertexPositions[i] << std::endl;
    	std::cout << "<" << pos(3*i) << ", " << pos(3*i + 1) << ", " << pos(3*i + 2) << ">" << std::endl;

    	assert(vpg.inputVertexPositions[i][0] == pos(3*i));
    	assert(vpg.inputVertexPositions[i][1] == pos(3*i+1));
    	assert(vpg.inputVertexPositions[i][2] == pos(3*i+2));
    }

    // Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>> positions(reinterpret_cast<double*>(d), 3, n_vertices);

	// Force f(*ptrmesh, *ptrvpg);

	// std::cout << "Sizeof Force: " << sizeof(Force) << std::endl;

	// //std::unique_ptr<force> ptrf(ptrmesh, ptrvpg);
	// Eigen::Matrix<double, Eigen::Dynamic, 3> bf = f.bending_force(1.0, 0);
	// std::cout << "bending force" << std::endl << bf << std::endl;
	// for (std::size_t i = 0; i < bf.rows(); ++i) {
	// 	std::cout << "force mag" <<  bf.row(i).norm() << std::endl;
	// }

	// gcs::IntrinsicGeometryInterface& geometry = *ptrvpg;

	// for (gcs::Face f : ptrmesh->faces()) {

	// 	// Managed array holding quantity
	// 	double area = geometry.faceAreas[f];

	// 	// Immediate computation, computes directly from
	// 	// input data without touching caches.
	// 	// Generally discouraged but occasionally useful.
	// 	area = ptrvpg->faceArea(f);
	// 	//std::cout << area << std::endl;
	// }

	// // Compute vertex area
	// gcs::VertexData<double> vertArea(*ptrmesh, 0.);
	// for (gcs::Vertex v : ptrmesh->vertices()) {
	// 	for (gcs::Face f : v.adjacentFaces()) {
	// 		vertArea[v] += geometry.faceAreas[f] / f.degree();
	// 	}
	// 	//std::cout << "degree =" << f.degree() << std::endl;
	// 	//std::cout << vertArea[v] << std::endl;
	// }


	// polyscope::init();
	// polyscope::registerSurfaceMesh("myMesh",
	// 							   ptrvpg->inputVertexPositions,
	//  							   ptrmesh->getFaceVertexList());
	// polyscope::show();

	Eigen::Matrix<double, 3, 3> foo;
	foo << 1, 2, 3,
  		   4, 5, 6,
           7, 8, 9;

    std::cout << foo << std::endl;

    foo.row(1) << 11, 12, 13;

    std::cout << foo << std::endl;

    auto tmp = foo.row(1);
    std::cout << "decltype(tmp) is " << type_name<decltype(tmp)>() << '\n';

    tmp << 14, 15, 16;
    std::cout << foo << std::endl;

	return 0;
}


