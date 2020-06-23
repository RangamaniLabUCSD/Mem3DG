#include "ddgsolver/force.h"
#include "ddgsolver/integrator.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <pcg_random.hpp>

#include <iostream>

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void velocityVerlet(Force &f, double dt, double total_time, double tolerance) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> force;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newForce;
  force.resize(f.mesh.nVertices(), 3);
  force.setZero();
  newForce.resize(f.mesh.nVertices(), 3);
  newForce.setZero();

  auto vel_e = ddgsolver::EigenMap<double, 3>(f.vel);
  auto pos_e = ddgsolver::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  const double hdt = 0.5 * dt;
  const double hdt2 = hdt * dt;

  for (int i = 0; i < total_time / dt; i++) {
    // Update all forces
    //f.getBendingForces();
    //f.getStretchingForces();
    //f.getPressureForces();
    f.getConservativeForces();
    f.getDPDForces();
    f.getExternalForces();

    pos_e +=
        (vel_e.rowwise() - (vel_e.colwise().sum() / f.mesh.nVertices())) * dt +
        force * hdt2;

    auto staticForce = ddgsolver::EigenMap<double, 3>(f.bendingForces) +
                       ddgsolver::EigenMap<double, 3>(f.stretchingForces) +
                       ddgsolver::EigenMap<double, 3>(f.pressureForces) +
                       ddgsolver::EigenMap<double, 3>(f.externalForces);
    auto dynamicForce = ddgsolver::EigenMap<double, 3>(f.dampingForces) +
                        ddgsolver::EigenMap<double, 3>(f.stochasticForces);
    auto newForce = staticForce + dynamicForce;

    vel_e += (force + newForce) * hdt;
    force = newForce;
    f.update_Vertex_positions(); // recompute cached values;

    double staticForce_mag = staticForce.norm();
    if (staticForce_mag < tolerance) {
      break;
    }
    // std::cout << "Force: " << staticForce_mag << std::endl;
    std::cout << "process: " << i << std::endl;
  }
}
} // namespace ddgsolver
