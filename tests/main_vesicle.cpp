#include <iostream>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;


int main() {
	/// geometric parameters
	int nSub = 3;

	/// physical parameters 
	ddgsolver::Parameters p;
	p.Kb = 0.01;			//Kb
	p.H0 = 1.4 * 0;				//H0
	p.Kse = 0;      //Kse
	p.Ksl = 1;				//Ksl
	p.Ksg = 2;				//Ksg
	p.Kv = 1;			  //Kv
	p.gamma = 1 * 0;				//gamma
	p.Vt = 0.7;			//Vt
	p.kt = 0.00001 * 0;		//Kt 

	p.ptInd = 0;       
	p.extF = 0.2 * 0;
	p.conc = 25 * 0;

	/// integration parameters
	double h = 0.001;
	double T = 200;
	double eps = 1e-9;// 1e-9;
	double tSave = 10; // save after time tSave

	p.sigma = sqrt(2 * p.gamma * p.kt / h);

	/// choose the starting mesh 
	std::string option = "sphere"; // 1. "sphere" 
																 // 2. "nameOfTheFile" = "input-file/%%%%%.ply"

	/// initialize mesh and vpg 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

	/// construct the starting mesh based on "option"
	if (option == "sphere") {
		/// initialize icosphere 
		std::vector<gc::Vector3> coords;
		std::vector<std::vector<std::size_t>> polygons;
		ddgsolver::icosphere(coords, polygons, nSub);
		gcs::PolygonSoupMesh soup(polygons, coords);
		soup.mergeIdenticalVertices();
		std::tie(ptrmesh, ptrvpg) =
		gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
	}
	else{
		std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(option);
	}

	auto& mesh = *ptrmesh;
	auto& vpg = *ptrvpg;

	gcs::RichSurfaceMeshData plyData(mesh);
	plyData.addGeometry(vpg);

	/// run the program based on "run"
	ddgsolver::Force f(mesh, vpg, p);
	ddgsolver::integrator integration(mesh, vpg, plyData, f, h, T, p, eps, tSave);
	integration.velocityVerlet();
	//velocityVerlet(f, h, T, eps, tSave);

	return 0;
	}
