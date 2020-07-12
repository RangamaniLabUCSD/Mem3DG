#include <iostream>

#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/utilities/vector3.h>

#include "ddgsolver/force.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;


int main() {

	std::string outputFolder = "output-file/";

	/// physical parameters 
	ddgsolver::Parameters p;
	p.Kb = 0.01;			//Kb
	p.H0 = 1.4 * 0;				//H0
	p.Kse = 0;      //Kse
	p.Ksl = 3;				//Ksl
	p.Ksg = 0;				//Ksg
	p.Kv = 1;			  //Kv
	p.gamma = 1;				//gamma
	p.Vt = 0.7;			//Vt
	p.kt = 0.00001;		//Kt 
	p.ptInd = 0;
	p.extF = 0.2;
	p.conc = 25;

	/// integration parameters
	double h = 0.0001;
	double T = 1;
	double eps = 1e-9;// 1e-9;
	double tSave = 0.25; // save after time tSave

	p.sigma = sqrt(2 * p.gamma * p.kt / h);

	/// choose the starting mesh 
	std::string option = "input-file/UVsphere.ply"; // 1. "input-file/UVsphere.ply" 
																								// 2. "input-file/Vt_%d_H0_%d.ply"

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
	std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(option);
	auto& mesh = *ptrmesh;
	auto& vpg = *ptrvpg;

	gcs::RichSurfaceMeshData plyData(mesh);
	plyData.addGeometry(vpg);

	/// run the program based on "run"
	ddgsolver::Force f(mesh, vpg, p);
	ddgsolver::integration::velocityVerlet(f, h, T, eps, tSave, outputFolder);

	return 0;
}
