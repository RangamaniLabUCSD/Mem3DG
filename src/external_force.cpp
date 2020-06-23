#include <iostream>

#include <Eigen/Core>

#include "ddgsolver/force.h"

namespace ddgsolver {

	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;

	void Force::getExternalForces() {
		auto externalForces_e = ddgsolver::EigenMap<double, 3>(externalForces);
		if (extF != 0) {
			// std::cout << "external force::" << appliedForceMagnitude << std::endl;
			// externalForce_e = vertexAreaGradientNormal.array().colwise() *
			//	appliedForceMagnitude.array();
			externalForces_e = appliedForceMagnitude * vertexAreaGradientNormal.row(ptInd);
		}
		else {
			externalForces_e.setZero();
		}
	}
}
