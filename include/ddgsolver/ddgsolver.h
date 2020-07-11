
#include <iostream>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/simple_polygon_mesh.h>

#include "ddgsolver/force.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/icosphere.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

int genIcosphere(size_t nSub, std::string path){
	/// initialize mesh and vpg 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;

	/// initialize icosphere 
	std::vector<gc::Vector3> coords;
	std::vector<std::vector<std::size_t>> polygons;
	ddgsolver::icosphere(coords, polygons, nSub);
	gcs::SimplePolygonMesh soup(polygons, coords);
	soup.mergeIdenticalVertices();
	std::tie(ptrmesh, ptrvpg) =
		gcs::makeHalfedgeAndGeometry(soup.polygons, soup.vertexCoordinates);
	writeSurfaceMesh(*ptrmesh, *ptrvpg, path);

	return 0;
}

int driver(std::string inputMesh, double Kb, double H0,
							double Kse, double Ksl, double Ksg,
							double Kv, double Vt, double gamma,
							double kt, size_t ptInd, double extF,
							double conc, double h, double T,
							double eps, double tSave) {

	/// physical parameters 
	double sigma = sqrt(2 * gamma * kt / h);
	ddgsolver::Parameters p{ Kb, H0,
		Ksl, Ksg, Kse,
		Kv, gamma, Vt,
	  kt, sigma, ptInd, extF,
		conc};

	//std::unique_ptr<gcs::SurfaceMesh> ptrmesh;
	//std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	//std::unique_ptr<gcs::RichSurfaceMeshData> richData;
	//std::tie(ptrmesh, richData) = gcs::RichSurfaceMeshData::readMeshAndData(option);
	//ptrvpg = richData->getGeometry();
	//auto& plyData = *richData;
	//auto& mesh = *ptrmesh;
	//auto& vpg = *ptrvpg;

	/// initialize mesh and vpg 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(inputMesh);
	auto& mesh = *ptrmesh;
	auto& vpg = *ptrvpg;

	gcs::RichSurfaceMeshData plyData(mesh);
	plyData.addGeometry(vpg);

	/// run the program based on "run"
	ddgsolver::Force f(mesh, vpg, p);
	ddgsolver::integration::velocityVerlet(f, h, T, eps, tSave);

	return 0;

}
