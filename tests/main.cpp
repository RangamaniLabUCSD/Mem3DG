
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// overload << to print vector;
template <typename T>
std::ostream &operator<<(std::ostream &output, const std::vector<T> &v) {
  output << "[";
  for (size_t i = 0; i != v.size() - 1; ++i) {
    output << v[i] << ",";
  }
  output << v[v.size() - 1];
  output << "]";
  return output;
}

int main() {

	//initialization of the geometry and mesh
	// Construct an icosahedron
	// Construct a std::vector of Vector 3
	std::vector<gc::Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons;

	ddgsolver::icosphere(coords, polygons, 2);

	gc::PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();

	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = 
		gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

	auto& mesh = *ptrmesh;
	auto& vpg  = *ptrvpg;
	
	//auto pos = mapVecToEigen<double, 3>(vpg.inputVertexPositions.rawdata());
	//std::cout << "decltype(pos) is " << type_name<decltype(pos)>() << '\n';
	/*for (std::size_t i = 0; i < mesh.nVertices(); ++i) {
		std::cout << vpg.inputVertexPositions[i] << std::endl;
		std::cout << "<" << pos(3 * i) << ", " << pos(3 * i + 1) << ", " << pos(3 * i + 2) << ">" << std::endl;
		assert(vpg.inputVertexPositions[i][0] == pos(3 * i));
		assert(vpg.inputVertexPositions[i][1] == pos(3 * i + 1));
		assert(vpg.inputVertexPositions[i][2] == pos(3 * i + 2));
	}*/
	

	struct parameters {
		double Kb, H0, Ksl, Ksg, Kv, gamma, Vt, kt;
	};
	
	parameters p{
							0.03,			//Kb
							0,				//H0
							4,				//Ksl
							6,				//Ksg
							1,				//Kv
							1,				//gamma
							0.6,			//Vt
							0.0001 };		//Kt     

	double h = 0.01;
	double T = 40;

	double sigma = sqrt(2 * p.gamma * p.kt / h);
	// initiate force object f
	ddgsolver::Force f(mesh,vpg,h);
	std::cout << "Sizeof Force: " << sizeof(f) << std::endl;

	gcs::FaceData<size_t> faceInd = vpg.faceIndices;
	gc::Vector3 totalForce;
	for (size_t i = 0; i < T / h; i++) {
		/*polyscope::registerSurfaceMesh("myMesh",
			ptrvpg->inputVertexPositions,
			ptrmesh->getFaceVertexList());*/
		//polyscope::show();
		f.getBendingForces(p.Kb, p.H0);
		f.getStretchingForces(p.Ksl, p.Ksg);
		f.getPressureForces(p.Kv, p.Vt);
		f.getDampingForces(p.gamma);
		f.getStochasticForces(sigma);

		//if (i == 18) { polyscope::show(); }
		//for (gcs::Vertex v : mesh.vertices()) {
		//	std::cout << "i: " << i << std::endl;
		//	std::cout << "v: " << v << std::endl;
		//	/*std::cout << "bending force" << i << f.bendingForces[v] << std::endl;
		//	std::cout << "stretching force" << f.stretchingForces[v] << std::endl;
		//	std::cout << "damping force" << f.dampingForces[v] << std::endl;
		//	std::cout << "pressure force" << f.pressureForces[v] << std::endl;
		//	std::cout << "stochastic force" << f.stochasticForces[v] << std::endl;*/
		//}
		//	

		gcs::VertexData<gc::Vector3> temp = vpg.inputVertexPositions;
		for (gcs::Vertex v : mesh.vertices()) {
			bool flag = true;
			for (gcs::Face f : v.adjacentFaces()) {
				if (faceInd[f] == 0) {
					flag = false;
				}
			}
			if (flag == true) {
				vpg.inputVertexPositions[v] *= 2;
				totalForce = f.bendingForces[v]
					+ f.stretchingForces[v]
					+ f.pressureForces[v]
					+ f.dampingForces[v]
					+ f.stochasticForces[v];
				vpg.inputVertexPositions[v] += totalForce * h * h - f.pastPositions[v];
			}	
		}
		std::cout << "totalForce:  " << totalForce.norm() << std::endl;
		f.update_Vertex_positions();
		f.pastPositions = temp;
	}
	polyscope::init();
	polyscope::registerSurfaceMesh("myMesh",
		ptrvpg->inputVertexPositions,
		ptrmesh->getFaceVertexList());
	polyscope::show();


	//f.getBendingForces(Kb, H0);
	//f.getStretchingForces(Ksl, Ksg);
	//f.getPressureForces(Kv, Vt);
	//f.getDampingForces(gamma);
	//f.getStochasticForces(sigma);
	//for (gcs::Vertex v : mesh.vertices()) {
	//	//std::cout << "bending force" << f.bendingForces[v] << std::endl;
	//	std::cout << "stretching force" << f.stretchingForces[v] << std::endl;
	//	//std::cout << "pressure force" << f.pressureForces[v] << std::endl;
	//	//std::cout << "damping force" << f.dampingForces[v] << std::endl;
	//	//std::cout << "stochastic force" << f.stochasticForces[v] << std::endl;
	//}


	// std::cout << "Sizeof Force: " << sizeof(Force) << std::endl;




	return 0;
}
