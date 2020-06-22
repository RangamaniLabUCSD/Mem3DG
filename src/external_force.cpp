#include <iostream>

#include <Eigen/Core>

#include "ddgsolver/force.h"

namespace ddgsolver {

	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;

	void Force::getExternalForces() {
		// alias distance
		auto externalForce_e = ddgsolver::EigenMap<double, 3>(externalForces);
		// std::cout << "external force::" << appliedForceMagnitude << std::endl;
		externalForce_e = vertexAreaGradientNormal.array().colwise() *
			appliedForceMagnitude.array();
	}
}
