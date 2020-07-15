#include <iostream>

#include <Eigen/Core>

#include "ddgsolver/force.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void Force::getExternalForces() {
  auto externalForces_e = EigenMap<double, 3>(externalForces);
  if (P.extF != 0) {
    // std::cout << "external force::" << appliedForceMagnitude << std::endl;
    // externalForce_e = vertexAreaGradientNormal.array().colwise() *
    //	appliedForceMagnitude.array();
    auto vertexAngleNormal_e = EigenMap<double, 3>(vpg.vertexNormals);
    externalForces_e = appliedForceMagnitude * vertexAngleNormal_e.row(P.ptInd);
  } else {
    externalForces_e.setZero();
  }
}
} // namespace ddgsolver
