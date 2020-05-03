#include "ddgsolver/force.h"
#include "geometrycentral/utilities/vector3.h"
#include <math.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include "geometrycentral/utilities/vector3.h"
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include "geometrycentral/numerical/linear_solvers.h"
//#define NDEBUG
#include <assert.h>  

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

double Force::signedVolumeFromFace(gcs::Face& f, gcs::VertexPositionGeometry& vpg) {
	
	gc::Vector3 p[3];
	size_t i = 0;
	for (gcs::Vertex v : f.adjacentVertices()) {
		p[i] = vpg.inputVertexPositions[v];
		i++;
	}
	double v321 = p[2].x * p[1].y * p[0].z;
	double v231 = p[1].x * p[2].y * p[0].z;
	double v312 = p[2].x * p[0].y * p[1].z;
	double v132 = p[0].x * p[2].y * p[1].z;
	double v213 = p[1].x * p[0].y * p[2].z;
	double v123 = p[0].x * p[1].y * p[2].z;

	return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123);
}


void Force::getPressureForces(double Kv, double Vt) {
	//gc::Vector3 origin { 0.0,0.0,0.0 };
	double total_volume = 0;
	double face_volume;
	gcs::FaceData<int> sign_of_volume(mesh);
	for (gcs::Face f : mesh.faces()) {
		face_volume = signedVolumeFromFace(f, vpg);
		total_volume += face_volume;
		if (face_volume < 0) {
			sign_of_volume[f] = -1;
		}
		else {
			sign_of_volume[f] = 1;
		}
	}
	//std::cout << "total volume" << total_volume << std::endl;

	for (gcs::Vertex v : mesh.vertices()) {
		for (gcs::Halfedge he : v.outgoingHalfedges()) {
			gcs::Halfedge base_he = he.next();
			gc::Vector3 p1 = vpg.inputVertexPositions[base_he.vertex()];
			gc::Vector3 p2 = vpg.inputVertexPositions[base_he.next().vertex()];
			gc::Vector3 dVdx = 0.5 * gc::cross(p1, p2) / 3.0;
			assert(gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) > 0);
			//std::cout << "i am here" << (gc::dot(dVdx, vpg.inputVertexPositions[v] - p1) < 0) <<  std::endl;
			dVdx *= sign_of_volume[he.face()];
			pressureForces[v] = -0.5 *Kv * (total_volume - targetVolume * Vt) / (targetVolume * Vt) * dVdx;
			//for (size_t i = 0; i < 3; i++) {
			//	pressureForces(v_ind[v], i) += dVdx[i];
			//}
			//force.row(v_ind[v]) << dVdx.x, dVdx.y, dVdx.z;
		}
	}
}