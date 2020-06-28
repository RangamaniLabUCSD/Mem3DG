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
  Eigen::Matrix<double, Eigen::Dynamic, 3> force;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newForce;
  force.resize(f.mesh.nVertices(), 3);
  force.setZero();
  newForce.resize(f.mesh.nVertices(), 3);
  newForce.setZero();

  int nSave = int(tSave / dt);

  auto vel_e = ddgsolver::EigenMap<double, 3>(f.vel);
  auto pos_e = ddgsolver::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  const double hdt = 0.5 * dt;
  const double hdt2 = hdt * dt;

  Eigen::Matrix<double, Eigen::Dynamic, 3> staticForce;
  Eigen::Matrix<double, Eigen::Dynamic, 3> dynamicForce;

  this -> getLogFiles();

  for (int i = 0; i < total_time / dt; i++) {
    // Update all forces
    //f.getBendingForces();
    //f.getStretchingForces();
    //f.getPressureForces();
    f.getConservativeForces();
    f.getDPDForces();
    f.getExternalForces();

    /*std::cout << "bf: " << ddgsolver::EigenMap<double, 3>(f.bendingForces).norm()
      << "sf: " << ddgsolver::EigenMap<double, 3>(f.stretchingForces).norm()
      << "pf: " << ddgsolver::EigenMap<double, 3>(f.pressureForces).norm()
      << "df: " << ddgsolver::EigenMap<double, 3>(f.dampingForces).norm()
      << "xf: " << ddgsolver::EigenMap<double, 3>(f.stochasticForces).norm() << std::endl;*/

    pos_e +=
      (vel_e.rowwise() - (vel_e.colwise().sum() / f.mesh.nVertices())) * dt +
      force * hdt2;

    staticForce = ddgsolver::EigenMap<double, 3>(f.bendingForces) +
      ddgsolver::EigenMap<double, 3>(f.stretchingForces) +
      ddgsolver::EigenMap<double, 3>(f.pressureForces) +
      ddgsolver::EigenMap<double, 3>(f.externalForces);
    dynamicForce = ddgsolver::EigenMap<double, 3>(f.dampingForces) +
      ddgsolver::EigenMap<double, 3>(f.stochasticForces);
    newForce = staticForce + dynamicForce;

    vel_e += (force + newForce) * hdt;
    force = newForce;
    f.update_Vertex_positions(); // recompute cached values;
    double staticForce_mag = staticForce.norm();

    if (staticForce_mag < tolerance) {
      break;
    }

    if ((i % nSave == 0) || (i == int(total_time / dt))) {
      gcs::PlyHalfedgeMeshData data(f.mesh);
      data.addGeometry(f.vpg);
      char buffer[50];
      sprintf(buffer, "output-file/t=%d.ply", int(i * dt * 100));
      data.write(buffer);
      std::cout << "time: " << i * dt << std::endl;
      std::cout << "force: " << staticForce_mag << std::endl;
      std::cout << "area: " << f.surfaceArea / f.initialSurfaceArea << std::endl;
      std::cout << "total volume:  " << f.volume / f.maxVolume / f.Vt << std::endl;
    }

  }

}

void velocityVerlet(Force &f, double dt, double total_time, double tolerance, double tSave) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> force;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newForce;
  force.resize(f.mesh.nVertices(), 3);
  force.setZero();
  newForce.resize(f.mesh.nVertices(), 3);
  newForce.setZero();

  int nSave = int(tSave / dt);

  auto vel_e = ddgsolver::EigenMap<double, 3>(f.vel);
  auto pos_e = ddgsolver::EigenMap<double, 3>(f.vpg.inputVertexPositions);

  const double hdt = 0.5 * dt;
  const double hdt2 = hdt * dt;

  Eigen::Matrix<double, Eigen::Dynamic, 3> staticForce;
  Eigen::Matrix<double, Eigen::Dynamic, 3> dynamicForce;

  //this -> getLogFiles();

  for (int i = 0; i < total_time / dt; i++) {
    // Update all forces
    //f.getBendingForces();
    //f.getStretchingForces();
    //f.getPressureForces();
    f.getConservativeForces();
    f.getDPDForces();
    f.getExternalForces();

    std::cout << "bf: " << ddgsolver::EigenMap<double, 3>(f.bendingForces).norm()
          << "sf: " << ddgsolver::EigenMap<double, 3>(f.stretchingForces).norm()
          << "pf: " << ddgsolver::EigenMap<double, 3>(f.pressureForces).norm()
          << "df: " << ddgsolver::EigenMap<double, 3>(f.dampingForces).norm()
          << "xf: " << ddgsolver::EigenMap<double, 3>(f.stochasticForces).norm() <<std::endl;

    pos_e +=
        (vel_e.rowwise() - (vel_e.colwise().sum() / f.mesh.nVertices())) * dt +
        force * hdt2;

    staticForce = ddgsolver::EigenMap<double, 3>(f.bendingForces) +
                       ddgsolver::EigenMap<double, 3>(f.stretchingForces) +
                       ddgsolver::EigenMap<double, 3>(f.pressureForces) +
                       ddgsolver::EigenMap<double, 3>(f.externalForces);
    dynamicForce = ddgsolver::EigenMap<double, 3>(f.dampingForces) +
                        ddgsolver::EigenMap<double, 3>(f.stochasticForces);
    newForce = staticForce + dynamicForce;

    vel_e += (force + newForce) * hdt;
    force = newForce;
    f.update_Vertex_positions(); // recompute cached values;

    double staticForce_mag = staticForce.norm();

    if (staticForce_mag < tolerance) {
      break;
    }

    if ((i % nSave == 0) || (i == int(total_time/dt))) {
      gcs::PlyHalfedgeMeshData data(f.mesh);
      data.addGeometry(f.vpg);
      char buffer[50];
      sprintf(buffer, "output-file/t=%d.ply", int(i * dt * 100));
      data.write(buffer);
      // std::cout << "Force: " << staticForce_mag << std::endl;
      std::cout << "time: " << i * dt << std::endl;
      std::cout << "force: " << staticForce_mag << std::endl;
      std::cout << "area: " << f.surfaceArea / f.initialSurfaceArea << std::endl;
      std::cout << "total volume:  " << f.volume / f.maxVolume / f.Vt << std::endl;
    }
  }
}
} // namespace ddgsolver
