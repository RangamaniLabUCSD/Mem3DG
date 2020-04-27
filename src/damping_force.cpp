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

void Force::damping_force(double gamma) {

	gcs::VertexData<size_t>& v_ind = vpg.vertexIndices;

	gcs::VertexData<gc::Vector3> velocity(mesh);
	for (gcs::Vertex v : mesh.vertices()) {
		velocity[v] = (vpg.inputVertexPositions[v] - vertex_position_past[v]) / time_step;
		//std::cout << "I am here" << vertex_position_past[v] << std::endl;
	}

	for (gcs::Vertex v : mesh.vertices()) {
		for (gcs::Vertex v_adj: v.adjacentVertices()) {
			gc::Vector3 velo_diff = velocity[v] - velocity[v_adj];
			gc::Vector3 posi_diff_unit = (vpg.inputVertexPositions[v] - vpg.inputVertexPositions[v_adj]).normalize();
			for (size_t i = 0; i < 3; i++) {
				df(v_ind[v], i) += (gc::dot(velo_diff, posi_diff_unit) * posi_diff_unit)[i];
			}
		}

	}
}