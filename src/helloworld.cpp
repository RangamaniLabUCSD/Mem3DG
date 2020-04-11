
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


using HalfedgeMesh = geometrycentral::surface::HalfedgeMesh;
using PolygonSoupMesh = geometrycentral::PolygonSoupMesh;
using Vector3 = geometrycentral::Vector3;


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

/*
Eigen::Matrix<double, Eigen::Dynamic, 1> Bending_force(Eigen::SparseMatrix<double> L, Eigen::SparseMatrix<double> M_inv){

}
*/


int main() {
	
//initialization of the geometry and mesh 
	// Construct an icosahedron 
	// Construct a std::vector of Vector 3
	std::vector<Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons; 

	icosphere(coords, polygons);
	
	PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();
	
	std::unique_ptr<HalfedgeMesh> mesh;
	std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg;
	std::tie(mesh, vpg) = geometrycentral::surface::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

	namespace geosurf = geometrycentral::surface;

	geosurf::IntrinsicGeometryInterface& geometry = *vpg;

	// try out some functionality of eigen
	std::cout << "Eigen version:" << EIGEN_MAJOR_VERSION << "." << std::endl;
	using namespace Eigen;
	Matrix3f A;
	Matrix3d B;
	
	Matrix<double, 10, 10> M1;
	A << 1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f;
	
	//std::cout << A << std::endl;
	std::cout << "B's size " << B.size() << std::endl;
	std::cout << "B's shape " << B.cols() << "x" << B.rows() << std::endl;
	
	for (int j = 0; j < B.rows(); ++j) {
		for (int i = 0; i < B.cols(); ++i) {
			B(i, j) = 0.0;
		}
	}
	std::cout << B << std::endl;

	// populate the quantity
	geometry.requireFaceAreas();
	geometry.requireCotanLaplacian();
	geometry.requireVertexGalerkinMassMatrix();
	geometry.requireVertexGaussianCurvatures();

	// weak(conformal) laplacian operator
	Eigen::SparseMatrix<double> L = geometry.cotanLaplacian;
	std::cout << Eigen::MatrixXd(L) << std::endl;

	// Mass matrix (Galerkin approximation)
	Eigen::SparseMatrix<double> M = geometry.vertexGalerkinMassMatrix;
	std::cout << Eigen::MatrixXd(L) << std::endl;

	// Gaussian curvature 
	geosurf::VertexData<double> KG = geometry.vertexGaussianCurvatures;
	Eigen::Matrix<double, Eigen::Dynamic, 1> v = KG.toVector();
	std::cout << "Gaussian" << v << std::endl; 


	for (geosurf::Face f : mesh->faces()) {

		// Managed array holding quantity
		double area = geometry.faceAreas[f];

		// Immediate computation, computes directly from 
		// input data without touching caches.
		// Generally discouraged but occasionally useful.
		area = vpg->faceArea(f);
		//std::cout << area << std::endl;
	}

	// Compute vertex area
	geosurf::VertexData<double> vertArea(*mesh, 0.);
	for (geosurf::Vertex v : mesh->vertices()) {
		for (geosurf::Face f : v.adjacentFaces()) {
			vertArea[v] += geometry.faceAreas[f] / f.degree();
		}
		//std::cout << "degree =" << f.degree() << std::endl;
		std::cout << vertArea[v] << std::endl;
	}


	// print properties of mesh
	size_t a;
	a = (mesh->nEdges());
	std::cout << a << std::endl;
	a = (mesh->nHalfedges());
	std::cout << a << std::endl;

	// visualization 
	/*polyscope::init();
	polyscope::registerSurfaceMesh("myMesh", vpg->inputVertexPositions,mesh->getFaceVertexList());
	polyscope::show();*/
	

	return 0;


}


