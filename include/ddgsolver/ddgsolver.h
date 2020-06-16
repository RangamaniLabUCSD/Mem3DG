
#include <iostream>

#include <geometrycentral/surface/halfedge_factories.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/ply_halfedge_mesh_data.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include "ddgsolver/force.h"
#include "ddgsolver/icosphere.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/typetraits.h"
#include "ddgsolver/util.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void driver(std::string run, std::string option, int nSub, double H0, double Vt,
	double h, double T, double eps) {
	/// physical parameters 
	ddgsolver::Parameters p;
	p.Kb = 0.01;			//Kb
	p.H0 = H0;				//H0
	p.Kse = 0.1;      //Kse
	p.Ksl = 1;				//Ksl
	p.Ksg = 2;				//Ksg
	p.Kv = 10;			  //Kv
	p.gamma = 1;				//gamma
	p.Vt = 1 * Vt;			//Vt
	p.kt = 0.00001;		//Kt 
}
