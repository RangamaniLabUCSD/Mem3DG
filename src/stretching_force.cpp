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

template<typename T>
void log(gcs::FaceData<T> face_a, gcs::HalfedgeMesh& mesh, std::string name) {
	for (gcs::Face f : mesh.faces()) {
		std::cout << name << face_a[f] << std::endl;
	}
}

gc::Vector3 vec_from_halfedge(gcs::Halfedge& he, gcs::VertexPositionGeometry& vpg) {
	gcs::Halfedge he_next = he.next();
	gc::Vector3 vec = vpg.inputVertexPositions[he_next.vertex()]
		- vpg.inputVertexPositions[he.vertex()];
	return vec;
}

Eigen::Matrix<double, Eigen::Dynamic, 3> force::stretching_force(double Ksl, double Ksg) {
	vpg.requireFaceNormals();
	gcs::FaceData<gc::Vector3>& face_n = vpg.faceNormals;
	log(face_n, mesh,"face normal");

	vpg.requireFaceAreas();
	gcs::FaceData<double>& face_a = vpg.faceAreas;
	log(face_a, mesh, "faceArea");

	vpg.requireVertexIndices();
	gcs::VertexData<size_t>& v_ind = vpg.vertexIndices;

	Eigen::Matrix<double, Eigen::Dynamic, 3> force;
	force.resize(mesh.nVertices(),3);

	for (gcs::Vertex v : mesh.vertices()) {
		
		for (gcs::Halfedge he : v.outgoingHalfedges()) {
			gcs::Halfedge base_he = he.next();
			gc::Vector3 base_vec = -vec_from_halfedge(base_he,vpg);
			//std::cout << "base vector" << base_vec << std::endl;
			gc::Vector3 gradient = gc::cross(base_vec, face_n[he.face()]);
			assert((gc::dot(gradient, vec_from_halfedge(he, vpg))) < 0);
			std::cout << "gradient" << gradient << std::endl;
			//auto force_v = force.row(v_ind[v]);
			//Eigen::Map<Eigen::Matrix<double, 1, 3>> force_v (&gradient.x, 3);
			force.row(v_ind[v]) << gradient.x, gradient.y, gradient.z;
		}
	}

	std::cout << "force" << force << std::endl;


	return force;
}