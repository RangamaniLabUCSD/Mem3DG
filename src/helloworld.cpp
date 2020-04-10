
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

	for (geosurf::Face f : mesh->faces()) {

		// Managed array holding quantity
		double area = geometry.faceAreas[f];

		// Immediate computation, computes directly from 
		// input data without touching caches.
		// Generally discouraged but occasionally useful.
		area = vpg->faceArea(f);
		std::cout << area << std::endl;
	}

	size_t a;
	a = (mesh->nEdges());
	std::cout << a << std::endl;
	a = (mesh->nHalfedges());
	std::cout << a << std::endl;

	//using namespace geometrycentral::surface
	polyscope::init();
	polyscope::registerSurfaceMesh("myMesh", vpg->inputVertexPositions,mesh->getFaceVertexList());
	polyscope::show();

	return 0;


}


//template<typename T>
//::std::ostream& operator<<(::std::ostream& output, const std::vector<T>& v) {
//    output << "[";
//    for (T item : v) {
//        output << item << ", ";
//    }
//    output << "]";
//    return output;
//}
//
//
//// strip unused vertices from face-vertex lists
//void stripUnusedVertices(std::vector<Vector3>& positions, std::vector<std::vector<size_t>>& faceIndices) {
//    size_t nVert = positions.size();
//
//
//    // Find any unused vertices
//    std::vector<size_t> vertexDegreeCount(nVert, 0);
//    size_t nUsedVerts = 0;
//    for (auto& f : faceIndices) {
//        for (auto& i : f) {
//            // Make sure we can safely index positions
//            GC_SAFETY_ASSERT(i < positions.size(),
//                "face index list has a vertex index which is greater than the number of vertices");
//            vertexDegreeCount[i]++;
//            if (vertexDegreeCount[i] == 1) {
//                nUsedVerts++;
//            }
//        }
//    }
//
//
//    // Early exit if dense
//    if (nUsedVerts == nVert) {
//        return;
//    }
//
//
//    // Else: strip unused vertices and re-index faces
//    size_t nNewVertices = 0;
//    std::vector<size_t> oldToNewVertexInd(nVert);
//    for (size_t iV = 0; iV < nVert; iV++) {
//        if (vertexDegreeCount[iV] > 0) {
//            oldToNewVertexInd[iV] = nNewVertices;
//            positions[nNewVertices] = positions[iV];
//            nNewVertices++;
//        }
//    }
//    positions.resize(nNewVertices);
//    for (auto& f : faceIndices) {
//        for (auto& i : f) {
//            i = oldToNewVertexInd[i];
//        }
//    }
//}
//
//
//
//int main() {
//    //std::unique_ptr<HalfedgeMesh> mesh = new HalfedgeMesh();
//    // Construct an icosahedron
//    // Construct a std::vector of Vector3
//    std::vector<Vector3> coords;
//    std::vector<std::vector<std::size_t>> polygons;
//    // [[1,2,3], [2,3,4],[..] ]
//    coords.push_back(Vector3{ 1.2f,0,0 });
//    coords.push_back(Vector3{ -1,0,0 });
//    coords.push_back(Vector3{ 0,1,0 });
//    coords.push_back(Vector3{ 0,0,1 });
//    for (Vector3 v : coords)
//        std::cout << v << std::endl;
//    polygons.push_back(std::vector<std::size_t>{0, 1, 2});
//    polygons.push_back(std::vector<std::size_t>{0, 1, 3});
//    polygons.push_back(std::vector<std::size_t>{0, 2, 3});
//    polygons.push_back(std::vector<std::size_t>{1, 2, 3});
//    for (std::vector<std::size_t> f : polygons)
//        std::cout << f << std::endl;
//    // std::cout << static_cast<std::size_t>(3.14159) << std::endl;
//    // int i = std::pow(2,31)-1;
//    // std::cout << i << ": " << std::bitset<32>(i) << std::endl;
//    // std::cout << i+1 << ": " << std::bitset<32>(i+1) << std::endl;
//    // std::size_t t = std::pow(2,31)-1;
//    // std::cout << t << ": " << std::bitset<32>(t) << std::endl;
//    // std::cout << t+1 << ": " << std::bitset<32>(t+1) << std::endl;
//    PolygonSoupMesh soup(polygons, coords);
//    soup.mergeIdenticalVertices();
//    // stripUnusedVertices(soup.vertexCoordinates, soup.polygons);
//    std::unique_ptr<HalfedgeMesh> mesh;
//    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> vpg;
//    std::tie(mesh, vpg) = geometrycentral::surface::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);
//    return 0;
//}