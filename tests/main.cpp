
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
	for (size_t i = 0; i < coords.size(); i++) {
		coords[i].x *= 1;
		//coords[i].y *= 3;
	}

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
	p.Kb = 0.01;			//Kb
	p.H0 = 3;				//H0
	p.Kse = 0.2;      //Kse
	p.Ksl = 1;				//Ksl
	p.Ksg = 2;				//Ksg
	p.Kv = 10;			//Kv
	p.gamma = 1;				//gamma
	p.Vt = 1 * 0.55;			//Vt
	p.kt = 0.00001;		//Kt 

	double h = 0.01;
	double T = 100;
	double eps = 1e-9;// 1e-9;

	
	ddgsolver::Force f(mesh,vpg);
	ddgsolver::integrator integration(mesh, vpg, f, h, T, p, eps);
	//integration.stormerVerlet();
	integration.velocityVerlet();

	//integration.p.H0 = 2;
	//integration.p.Vt = 0.8;
	//integration.velocityVerlet();

	polyscope::init();
	polyscope::registerSurfaceMesh("myMesh",
		ptrvpg->inputVertexPositions,
		ptrmesh->getFaceVertexList());
	polyscope::show();

	return 0;
	}

