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
#include "polyscope/curve_network.h"

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"
#include "ddgsolver/integrator.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

// parsing the option for two type of visualization
int main(int argc, char* argv[]) {

	if (argc < 3) {
		// Tell the user how to run the program
		std::cerr << "Usage: " << argv[0] << 
			" <typeOfMesh: 1. surface 2. network >" << " <.ply file>" << std::endl;
		return 1;
	}

	/// parameter needed for constructing the force object 
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
	p.extF = 0.2 * 0;
	p.conc = 25;
	double h = 0.0001;
	p.sigma = sqrt(2 * p.gamma * p.kt / h);

	/// type of mesh
	std::string typeMesh = argv[1];

	/// choose the .ply file 
	std::string option = argv[2]; // 1. "input-file/sphere.ply" 2. "input-file/Vt_%d_H0_%d.ply"

	/// initialize mesh and vpg 
	std::unique_ptr<gcs::HalfedgeMesh> ptrmesh;
	std::unique_ptr<gcs::VertexPositionGeometry> ptrvpg;
	std::unique_ptr<gcs::PlyHalfedgeMeshData> ptrPlyData;
	std::tie(ptrmesh, ptrvpg) = gcs::loadMesh(option);
	auto& mesh = *ptrmesh;
	auto& vpg = *ptrvpg;

	gcs::VertexData<double> meanCurvature = ptrPlyData->getVertexProperty<double>("mean curvature");
	gcs::VertexData<double> extForce = ptrPlyData->getVertexProperty<double>("external force");
	gcs::VertexData<double> normalForce = ptrPlyData->getVertexProperty<double>("fn");
	gcs::VertexData<double> tangentialForce = ptrPlyData->getVertexProperty<double>("ft");

	Eigen::Matrix<double, Eigen::Dynamic, 1> meanCurvature_e = meanCurvature.toMappedVector();
	Eigen::Matrix<double, Eigen::Dynamic, 1> extForce_e = extForce.toMappedVector();
	Eigen::Matrix<double, Eigen::Dynamic, 1> normalForce_e = normalForce.toMappedVector();
	Eigen::Matrix<double, Eigen::Dynamic, 1> tangentialForce_e = tangentialForce.toMappedVector();

	/// initialize force object f
	ddgsolver::Force f(mesh, vpg, p);
	f.getConservativeForces();

	/// Visualization: surface plot OR wire plot 
	polyscope::init();

	if (typeMesh == "surface") {

		polyscope::registerSurfaceMesh("Vesicle surface",
			ptrvpg->inputVertexPositions,
			ptrmesh->getFaceVertexList());
		polyscope::getSurfaceMesh("Vesicle surface")->addVertexScalarQuantity("mean curvature", meanCurvature_e);
		polyscope::getSurfaceMesh("Vesicle surface")->addVertexScalarQuantity("applied force", extForce_e);
		polyscope::getSurfaceMesh("Vesicle surface")->addVertexScalarQuantity("tangential force", tangentialForce_e);
		polyscope::getSurfaceMesh("Vesicle surface")->addVertexScalarQuantity("normal force", normalForce_e);
		polyscope::show();
	}

	else if (typeMesh == "network") {

		polyscope::registerCurveNetwork("Vesicle network",
			ptrvpg->inputVertexPositions,
			ptrmesh->getFaceVertexList());
		polyscope::getCurveNetwork("Vesicle surface")->addNodeScalarQuantity("mean curvature", meanCurvature_e);
		polyscope::getCurveNetwork("Vesicle surface")->addNodeScalarQuantity("applied force", extForce_e);
		polyscope::getCurveNetwork("Vesicle surface")->addNodeScalarQuantity("tangential force", tangentialForce_e);
		polyscope::getCurveNetwork("Vesicle surface")->addNodeScalarQuantity("normal force", normalForce_e);
		polyscope::show();
	}

	/// print message on polyscope and (screenshot)
	//char buffer[50];
	//sprintf(buffer, "Vt = %.2f, H0 = %.2f", p.Vt, p.H0);
	//polyscope::info(buffer);
	//sprintf(buffer, "output-file/Vt_%d_H0_%d.png", int(p.Vt * 100), int(p.H0 * 100));
	//polyscope::screenshot(buffer, true);

	return 0;
	}
