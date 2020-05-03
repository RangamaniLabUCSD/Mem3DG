
#include <iostream>
#include <math.h>
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

	icosphere(coords, polygons, 1);

	gc::PolygonSoupMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();

	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);

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
	
	double Kb = 0.01;
	double H0 = 0;
	double Ksl = 3;
	double Ksg = 6;
	double Kv = 1;
	double gamma = 1;
	double Vt = 0.7;
	double kt = 0.0001;

	double h = 0.01;
	double T = 1;

	double sigma = sqrt(2 * gamma * kt / h);
	// initiate force object f
	Force f(*ptrmesh, *ptrvpg,h);
	std::cout << "Sizeof Force: " << sizeof(f) << std::endl;

	polyscope::init();
	for (size_t i = 0; i < T / h; i++) {
		polyscope::registerSurfaceMesh("myMesh",
			ptrvpg->inputVertexPositions,
			ptrmesh->getFaceVertexList());
		//polyscope::show();
		f.getBendingForces(Kb, H0);
		f.getStretchingForces(Ksl, Ksg);
		f.getPressureForces(Kv, Vt);
		f.getDampingForces(gamma);
		f.getStochasticForces(sigma);

		if (i == 33) { polyscope::show(); }
		for (gcs::Vertex v : mesh.vertices()) {
			std::cout << "i: " << i << std::endl;
			std::cout << "v: " << v << std::endl;
			/*std::cout << "bending force" << i << f.bendingForces[v] << std::endl;
			std::cout << "stretching force" << f.stretchingForces[v] << std::endl;
			std::cout << "damping force" << f.dampingForces[v] << std::endl;
			std::cout << "pressure force" << f.pressureForces[v] << std::endl;
			std::cout << "stochastic force" << f.stochasticForces[v] << std::endl;*/
		}
			

		gcs::VertexData<Vector3> temp = vpg.inputVertexPositions;
		for (gcs::Vertex v : mesh.vertices()) {
			vpg.inputVertexPositions[v] *= 2;
			vpg.inputVertexPositions[v] += (f.bendingForces[v]
				+ f.stretchingForces[v]
				+ f.pressureForces[v]
				+ f.dampingForces[v]
				+ f.stochasticForces[v]) * h * h - f.pastPositions[v];
		}
		f.pastPositions = temp;
	}


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


