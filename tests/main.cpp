
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/ply_halfedge_mesh_data.h>
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
	/// physical parameters 
	ddgsolver::parameters p;
	p.Kb = 0.01;			//Kb
	p.H0 = 1.5;				//H0
	p.Kse = 0.1;      //Kse
	p.Ksl = 1;				//Ksl
	p.Ksg = 2;				//Ksg
	p.Kv = 10;			//Kv
	p.gamma = 1;				//gamma
	p.Vt = 1 * 0.7;			//Vt
	p.kt = 0.00001;		//Kt 

	// choose the run
	std::string run = "visualization"; // 1. "integration 2. "visualization

	/// Choose the starting mesh 
	std::string option = "output-file/Vt_70_H0_0.ply"; // 1. "sphere" 2. "continue" 3. "nameOfTheFile" = "output-file/Vt_%d_H0_%d.ply"

	/// initialize mesh and vpg 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

	if (option == "continue") {
		char buffer[50];
		sprintf(buffer, "output-file/Vt_%d_H0_%d.ply", int(p.Vt * 100), int(p.H0 * 100));
		std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(buffer);
	}
	else if (option == "sphere") {
		/// initialize icosphere 
		std::vector<gc::Vector3> coords;
		std::vector<std::vector<std::size_t>> polygons;
		ddgsolver::icosphere(coords, polygons, 2);
		gc::PolygonSoupMesh soup(polygons, coords);
		soup.mergeIdenticalVertices();
		std::tie(ptrmesh, ptrvpg) =
		gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates, true);
	}
	else {
		std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(option);
	}
	
	auto& mesh = *ptrmesh;
	auto& vpg = *ptrvpg;

	if (run == "integration") {
		/// integration parameters
		double h = 0.005;
		double T = 300;
		double eps = 1e-9;// 1e-9;

		/// solve
		ddgsolver::Force f(mesh, vpg);
		ddgsolver::integrator integration(mesh, vpg, f, h, T, p, eps);
		//integration.stormerVerlet();
		integration.velocityVerlet();

		/// save the .ply file  
		gcs::PlyHalfedgeMeshData data(mesh);
		data.addGeometry(vpg);
		char buffer[50];
		sprintf(buffer, "output-file/Vt_%d_H0_%d.ply", int(p.Vt * 100), int(p.H0 * 100));
		data.write(buffer);
	}


	/// visualization 
	polyscope::init();
	polyscope::registerSurfaceMesh("myMesh",
		ptrvpg->inputVertexPositions,
		ptrmesh->getFaceVertexList());
	polyscope::show();

	return 0;
	}

