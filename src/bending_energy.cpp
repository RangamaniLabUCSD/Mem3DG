#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>


namespace ddgsolver {
	namespace integration {

		namespace gc = ::geometrycentral;
		namespace gcs = ::geometrycentral::surface;

		double getTotalEnergy(Force& f) {
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

			/// bending energy 
			Eigen::Matrix<double, Eigen::Dynamic, 1> k_dH_sqrd;
			k_dH_sqrd.resize(f.Hn.rows(), 1);
			auto difference = f.Hn - f.H0n;
			k_dH_sqrd = f.P.Kb * (difference.array() * difference.array()).colwise().sum();
			double bE = (f.M * k_dH_sqrd).sum();

			/// stretching energy 
			//double sE =

			return bE;
		}
	} // namespace integration
} // namespace ddgsolver
