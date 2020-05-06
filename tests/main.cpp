
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
	
	ddgsolver::parameters p;
	p.Kb = 0.02;			//Kb
	p.H0 = 0;				//H0
	p.Ksl = 4;				//Ksl
	p.Ksg = 6;				//Ksg
	p.Kv = 1;			//Kv
	p.gamma = 1;				//gamma
	p.Vt = 0.6;			//Vt
	p.kt = 0.0001;		//Kt     

	double h = 0.01;
	double T = 20;

	p.sigma = sqrt(2 * p.gamma * p.kt / h);
	ddgsolver::Force f(mesh,vpg);
	ddgsolver::integrator integration(mesh, vpg, f, h, T, p);
	//std::cout << "hello I am here" << std::endl;
	//integration.stormerVerlet();
	integration.velocityVerlet();

	polyscope::init();
	polyscope::registerSurfaceMesh("myMesh",
		ptrvpg->inputVertexPositions,
		ptrmesh->getFaceVertexList());
	polyscope::show();

	return 0;
	}

