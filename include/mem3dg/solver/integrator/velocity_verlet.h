// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2021:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#pragma once

#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {

/**
 * @brief Velocity Verlet time Integration
 */
class DLL_PUBLIC VelocityVerlet : public Integrator {
private:
  // previous force
  gcs::VertexData<gc::Vector3> pastMechanicalForceVec;
  // total energy of the system
  double initialTotalEnergy;

public:
  bool isCapEnergy = true;
  VelocityVerlet(System &system_, double characteristicTimeStep_,
                 double totalTime_, double savePeriod_, double tolerance_,
                 std::string outputDirectory_)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_) {

    // print to console
    std::cout << "Running Velocity Verlet integrator ..." << std::endl;

    // check the validity of parameter
    checkParameters();

    pastMechanicalForceVec =
        gc::VertexData<gc::Vector3>(*system.mesh, {0, 0, 0});

    initialTotalEnergy = system.computeTotalEnergy();
    initialTotalEnergy -= system.energy.proteinInteriorPenalty;
  }

  /**
   * @brief velocity Verlet driver function
   */
  bool integrate() override;
  /**
   * @brief velocity Verlet marcher
   */
  void march() override;

  /**
   * @brief Velocity Verlet status computation and thresholding
   */
  void status() override;

  /**
   * @brief Check parameters for time integration
   */
  void checkParameters() override;

  /**
   * @brief step for n iterations
   */
  void step(std::size_t n) {
    for (std::size_t i = 0; i < n; i++) {
      status();
      march();
    }
  }
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
