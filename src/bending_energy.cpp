#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"
#include "ddgsolver/util.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>


namespace ddgsolver {
	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;

	void integrator::getTotalEnergy() {
		// comment: this may not be useful, the convergence can be be tested by checking its derivative 
		// which is the forces excluding the DPD forces. The energy trajectory of could actually numerically
		// integrated by post processing after saving all forces during the iterations. 

		//auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);
		//auto A = f.L.transpose() * f.M_inv * f.L; //could further cached since the same for all time pt
		//for (size_t i = 0; i < positions.rows(); i++) {
		//	for (size_t j = 0; j < positions.rows(); j++) {
		//		Eb += positions.row(i) * A * positions.transpose().col(j);
		//	}
		//}

		pastTotalEnergy = totalEnergy;

		/// bending energy 
		auto H_e = EigenMap(f.H);
		auto H0_e = EigenMap(f.H0);
		auto A = EigenMap(vpg.vertexDualAreas);
		Eigen::Matrix<double, Eigen::Dynamic, 1> difference = H_e - H0_e;
		Eigen::Matrix<double, Eigen::Dynamic, 1> k_dH_sqrd 
			= p.Kb * (difference.array() * difference.array());
		double bE = (A * k_dH_sqrd).sum();
		totalEnergy = bE;
		/// stretching energy 
		//double sE = 
	}
}

