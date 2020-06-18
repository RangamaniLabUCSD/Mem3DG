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

void integrator::velocityVerlet() {
  f.vel.fill({0, 0, 0});
  gcs::VertexData<gc::Vector3> totalForce(mesh, {0, 0, 0});
  gcs::VertexData<gc::Vector3> newTotalForce(mesh, {0, 0, 0});
  gc::Vector3 COMVelocity;
  Eigen::Map<Eigen::Matrix<double, 1, 3>> COMVelocity_e(&COMVelocity[0]);
  // auto COMVelocity_e = ddgsolver::EigenMap<double, 3>(COMVelocity);

  std::default_random_engine random_generator;
  std::uniform_int_distribution<int> uniform_dist(1, 10);

  for (size_t i = 0; i < timeSpan / timeStep; i++) {

    double timeStepHere = timeStep;
    // timeStep * (1 / pow(i, 1 / 2)); //* (1 / double(i)); //pow(i, 1 / 2));
    p.sigma = sqrt(2 * p.gamma * p.kt / timeStepHere);

    int prob = 11;

    int flag = uniform_dist(random_generator);
    if (true) {
      f.getBendingForces();
    }

    flag = uniform_dist(random_generator);
    if (flag < prob) {
      f.getStretchingForces();
    } else {
      f.stretchingForces.fill({0.0, 0.0, 0.0});
    }

    flag = uniform_dist(random_generator);
    if (flag < prob) {
      f.getPressureForces();
    } else {
      f.pressureForces.fill({0.0, 0.0, 0.0});
    }

    flag = uniform_dist(random_generator);
    if (flag < prob) {
      f.getDampingForces();
      f.getStochasticForces();
    } else {
      f.dampingForces.fill({0.0, 0.0, 0.0});
      f.stochasticForces.fill({0.0, 0.0, 0.0});
    }
    // f.getBendingForces(p.Kb, p.H0);
    // f.getStretchingForces(p.Ksl, p.Ksg,p.Kse);
    // f.getPressureForces(p.Kv, p.Vt);
    // f.getDampingForces(p.gamma);
    // f.getStochasticForces(p.sigma);

    COMVelocity_e =
        ddgsolver::EigenMap<double, 3>(f.vel).colwise().sum() /
        mesh.nVertices();
    for (gcs::Vertex v : mesh.vertices()) {
      vpg.inputVertexPositions[v] +=
          (f.vel[v] - COMVelocity) * timeStepHere +
          totalForce[v] * timeStepHere * timeStepHere * 0.5;

      newTotalForce[v] = f.bendingForces[v] + f.stretchingForces[v] +
                         f.pressureForces[v] + f.dampingForces[v] +
                         f.stochasticForces[v];

      /*std::cout << "bf: " << f.bendingForces[v].norm()
              << "sf: " << f.stretchingForces[v].norm()
              << "pf: " << f.pressureForces[v].norm()
              << "df: " << f.dampingForces[v].norm()
              << "xf: " << f.stochasticForces[v].norm() <<std::endl;*/

      f.vel[v] +=
          (totalForce[v] + newTotalForce[v]) * timeStepHere * 0.5;
      totalForce[v] = newTotalForce[v];
    }
    f.update_Vertex_positions();

    getBendingEnergy();
    if (((abs(pastBendingEnergy - bendingEnergy) / bendingEnergy) <
         tolerance) &&
        (i > 1) &&
        (abs(f.volume - f.targetVolume * p.Vt) / (f.targetVolume * p.Vt) <
         1e-2) &&
        (abs(f.surfaceArea - f.targetSurfaceArea) / (f.targetSurfaceArea) <
         1e-2)) {
      break;
    }
    // std::cout << "energy: " << bendingEnergy << std::endl;

    // std::cout << "process: " << int(double(i) / (timeSpan / timeStep) * 100)
    // << "%" << std::endl;
    std::cout << "process: " << i << std::endl;
  }
}

void velocityVerlet(Force &f, double dt, double total_time, double tolerance) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> force;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newForce;
  force.resize(f.mesh.nVertices(), 3);
  force.setZero();
  newForce.resize(f.mesh.nVertices(), 3);
  newForce.setZero();

  auto vel_e = ddgsolver::EigenMap<double, 3>(f.vel);
  auto pos_e = ddgsolver::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  f.sigma = sqrt(2 * f.gamma * f.kt / dt);
  const double hdt = 0.5 * dt;
  const double hdt2 = hdt * dt;

  for (int i = 0; i < total_time / dt; i++) {
    // Update all forces
    f.getBendingForces();
    f.getStretchingForces();
    f.getPressureForces();
    f.getDPDForces();
    // f.getDampingForces();
    // f.getStochasticForces();

    pos_e += (vel_e.rowwise() -
              (vel_e.colwise().sum() / f.mesh.nVertices())) *
                 dt +
             force * hdt2;

    newForce = EigenMap<double, 3>(f.bendingForces) +
               EigenMap<double, 3>(f.stretchingForces) +
               EigenMap<double, 3>(f.pressureForces) +
               EigenMap<double, 3>(f.dampingForces) +
               EigenMap<double, 3>(f.stochasticForces);

    vel_e += (force + newForce) * hdt;
    force.swap(newForce);
    f.update_Vertex_positions(); // recompute cached values;

    // getBendingEnergy();
    // if (((abs(pastBendingEnergy - bendingEnergy) / bendingEnergy) <
    //      tolerance) &&
    //     (i > 1) &&
    //     (abs(f.volume - f.targetVolume * f.Vt) / (f.targetVolume * f.Vt) <
    //      1e-2) &&
    //     (abs(f.surfaceArea - f.targetSurfaceArea) / (f.targetSurfaceArea) <
    //      1e-2)) {
    //   break;
    // }
    // std::cout << "energy: " << bendingEnergy << std::endl;

    std::cout << "BendingForces: " << dot(EigenMap<double, 3>(f.bendingForces), EigenMap<double, 3>(f.bendingForces)).cwiseSqrt().sum() << std::endl;
    std::cout << "StretchingForces: " << dot(EigenMap<double, 3>(f.stretchingForces), EigenMap<double, 3>(f.stretchingForces)).cwiseSqrt().sum() << std::endl;
    std::cout << "pressureForces: " << dot(EigenMap<double, 3>(f.pressureForces), EigenMap<double, 3>(f.pressureForces)).cwiseSqrt().sum() << std::endl;
    std::cout << "KE: " << 0.5*dot(vel_e, vel_e).cwiseSqrt().sum() << std::endl;
    std::cout << "process: " << i << std::endl;
  }
}
} // namespace ddgsolver
