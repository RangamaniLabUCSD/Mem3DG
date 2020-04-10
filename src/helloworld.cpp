
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

int main() {
	
//initialization of the geometry and mesh 
	// Construct an icosahedron 
	// Construct a std::vector of Vector 3
	std::vector<Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons; 

	// initialize the vertex coordinate
	coords.push_back(Vector3{ 0.0f, 0.0f, 2.0f });
	coords.push_back(Vector3{ 1.788854f, 0.000000f, 0.894427f });
	coords.push_back(Vector3{ 0.552786f, 1.701302f, 0.894427f });
	coords.push_back(Vector3{ -1.447214f, 1.051462f, 0.894427f });
	coords.push_back(Vector3{ -1.447214f, -1.051462f, 0.894427f });
	coords.push_back(Vector3{ 0.552786f, - 1.701302f, 0.894427f });
	coords.push_back(Vector3{ 1.447214f, 1.051462f, - 0.894427f});
	coords.push_back(Vector3{ -0.552786f, 1.701302f, - 0.894427f });
	coords.push_back(Vector3{ -1.788854f, 0.000000f, - 0.894427f });
	coords.push_back(Vector3{ -0.552786f, - 1.701302f, - 0.894427f });
	coords.push_back(Vector3{ 1.447214f, -1.051462f, - 0.894427f });
	coords.push_back(Vector3{ 0.0f , 0.0f,  - 2.0f });
	for (Vector3 v : coords) {
		std::cout << v << std::endl;
	}

	// initialize the face
	polygons.push_back(std::vector<std::size_t>{2, 0, 1});
	polygons.push_back(std::vector<std::size_t>{3, 0, 2	});
	polygons.push_back(std::vector<std::size_t>{4, 0, 3});
	polygons.push_back(std::vector<std::size_t>{5, 0, 4	});

	polygons.push_back(std::vector<std::size_t>{1, 0, 5});
	polygons.push_back(std::vector<std::size_t>{2, 1, 6	});
	polygons.push_back(std::vector<std::size_t>{7, 2, 6});
	polygons.push_back(std::vector<std::size_t>{3, 2, 7	});
		
	polygons.push_back(std::vector<std::size_t>{8, 3, 7});
	polygons.push_back(std::vector<std::size_t>{4, 3, 8	});
	polygons.push_back(std::vector<std::size_t>{9, 4, 8});
	polygons.push_back(std::vector<std::size_t>{5, 4, 9	});

	polygons.push_back(std::vector<std::size_t>{10, 5, 9});
	polygons.push_back(std::vector<std::size_t>{6, 1, 10	});
	polygons.push_back(std::vector<std::size_t>{1, 5, 10});
	polygons.push_back(std::vector<std::size_t>{6, 11, 7	});

	polygons.push_back(std::vector<std::size_t>{7, 11, 8});
	polygons.push_back(std::vector<std::size_t>{8, 11, 9	});
	polygons.push_back(std::vector<std::size_t>{9, 11, 10});
	polygons.push_back(std::vector<std::size_t>{10, 11, 6	});
	for (std::vector<std::size_t> f: polygons) {
		std::cout << f << std::endl;
	}
	
	PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();
	std::unique_ptr<HalfedgeMesh> mesh;
	std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg;
	std::tie(mesh, vpg) = geometrycentral::surface::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);


	namespace geosurf = geometrycentral::surface;
	geosurf::IntrinsicGeometryInterface& geometry = *vpg;

	// populate the quantity
	geometry.requireFaceAreas();
	geometry.requireCotanLaplacian();
	geometry.requireVertexGalerkinMassMatrix();
	geometry.requireVertexGaussianCurvatures();

	Eigen::SparseMatrix<double> L = geometry.cotanLaplacian;
	std::cout << Eigen::MatrixXd(L) << std::endl;

	L = geometry.vertexGalerkinMassMatrix;
	std::cout << Eigen::MatrixXd(L) << std::endl;

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


	size_t a;
	a = (mesh->nEdges());
	std::cout << a << std::endl;
	a = (mesh->nHalfedges());
	std::cout << a << std::endl;

	/*polyscope::init();
	polyscope::registerSurfaceMesh("myMesh", vpg->inputVertexPositions,mesh->getFaceVertexList());
	polyscope::show();*/
	

	return 0;


}


