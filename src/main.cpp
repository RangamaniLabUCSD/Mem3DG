
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

	//initialization of the geometry and mesh
	// Construct an icosahedron
	// Construct a std::vector of Vector 3
	std::vector<gc::Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons;

	icosphere(coords, polygons, 1);

	gc::PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();

	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

	force f(*ptrmesh, *ptrvpg);

	std::cout << "Sizeof Force: " << sizeof(force) << std::endl;

	//std::unique_ptr<force> ptrf(ptrmesh, ptrvpg);
	Eigen::Matrix<double, Eigen::Dynamic, 3> bf = f.bending_force(1.0, 0);
	std::cout << "bending force" << std::endl << bf << std::endl;
	std::cout << "no of vertices" << bf.rows() << std::endl;
	for (std::size_t i = 0; i < bf.rows(); ++i) {
		std::cout << "force mag" <<  bf.row(i).norm() << std::endl;
	}

	gcs::IntrinsicGeometryInterface& geometry = *ptrvpg;

	for (gcs::Face f : ptrmesh->faces()) {

		// Managed array holding quantity
		double area = geometry.faceAreas[f];

		// Immediate computation, computes directly from
		// input data without touching caches.
		// Generally discouraged but occasionally useful.
		area = ptrvpg->faceArea(f);
		//std::cout << area << std::endl;
	}

	// Compute vertex area
	gcs::VertexData<double> vertArea(*ptrmesh, 0.);
	for (gcs::Vertex v : ptrmesh->vertices()) {
		for (gcs::Face f : v.adjacentFaces()) {
			vertArea[v] += geometry.faceAreas[f] / f.degree();
		}
		//std::cout << "degree =" << f.degree() << std::endl;
		//std::cout << vertArea[v] << std::endl;
	}

	/*
	// print properties of ptrmesh
	size_t a;
	a = (ptrmesh->nEdges());
	std::cout << a << std::endl;
	a = (ptrmesh->nHalfedges());
	std::cout << a << std::endl;
	*/


	/// Explicit COPY based approach
	// gcs::VertexData<Vector3>& ivp = ptrvpg->inputVertexPositions;
	// Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1> eivp = ptrvpg->inputVertexPositions.toVector();

	// //manipulate eivp
	// // ...
	// //then copy it back
	// ptrvpg->fromVector(eivp);

	// std::cout << "EIVP: " << eivp << std::endl;


	/// Compiler error no access to protected member object.
	// auto & bar = ptrvpg->inputVertexPositions.data;


	/// CAN WE AVOID THE COPY?
	// .rawdata() is the function that we added to access the protected member "data" 
	// vector of vector3 
	std::vector<gc::Vector3, Eigen::aligned_allocator<gc::Vector3>>& foo = ptrvpg->inputVertexPositions.rawdata();
	for (int i = 0; i < foo.size(); ++i) {
		std::cout << "Foooo[" << i << "]: " << foo[i] << std::endl;

	}
	
	// [x,y,z;x,y,z;x,y,z;]

	// .data() returns a pointer to the first element in the array used internally by the vector.
	gc::Vector3* d = ptrvpg->inputVertexPositions.rawdata().data();

	Eigen::Map<Eigen::Matrix<gc::Vector3, Eigen::Dynamic, 1>> evec3(d, foo.size());
	

	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1 >> evecdouble(reinterpret_cast<double*>(d), foo.size()*3);
	
	for (int i = 0; i < foo.size(); ++i){
		std::cout << "Foo[" << i << "]: " << foo[i] << std::endl;


		std::cout << "evecdouble: " << evecdouble(3*i) << ", " << evecdouble(3*i+1) << ", " << evecdouble(3*i+2) << std::endl;
	}
	 //std::cout << evec3 << std::endl;
	 //std::cout << evecdouble << std::endl;

	 int i = 0;
	 for(auto& v : foo){
	 	 std::cout << v << std::endl;
	 	v[0] *= 5;
	 	v[1] *= 1;
	 	v[2] *= 1;

	 	std::cout << "bfr: " << v << std::endl;
	 	std::cout << "cmp: " << ptrvpg->inputVertexPositions[i] << std::endl;
	 	++i;
	 }

	 /*polyscope::init();
	 polyscope::registerSurfacemesh("myMesh",
	 							   ptrvpg->inputVertexPositions,
	 							   ptrmesh->getFaceVertexList());
	 polyscope::show();*/

	return 0;
}


