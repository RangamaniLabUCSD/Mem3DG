#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>


namespace ddgsolver {
	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;

	void integrator::getBendingEnergy() {
		//auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);
		//auto A = f.L.transpose() * f.M_inv * f.L; //could further cached since the same for all time pt
		//for (size_t i = 0; i < positions.rows(); i++) {
		//	for (size_t j = 0; j < positions.rows(); j++) {
		//		Eb += positions.row(i) * A * positions.transpose().col(j);
		//	}
		//}
		pastBendingEnergy = bendingEnergy;

		Eigen::Matrix<double, Eigen::Dynamic, 1> k_dH_sqrd;
		k_dH_sqrd.resize(f.Hn.rows(), 1);

		auto difference = f.Hn - f.H0n;

		for (size_t row = 0; row < f.Hn.rows(); row++) {
			k_dH_sqrd(row) = p.Kb * (difference.row(row).dot(difference.row(row)));
		}

		bendingEnergy = (f.M * k_dH_sqrd).sum();
	}
}

