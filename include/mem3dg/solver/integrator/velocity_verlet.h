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
  // total pressure
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure;
  // total pressure of new iteration
  Eigen::Matrix<double, Eigen::Dynamic, 3> newTotalPressure;
  // total energy of the system
  double totalEnergy;

public:
  VelocityVerlet(System &f_, double dt_, double total_time_, double tSave_,
                 double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_) {

    // print to console
    if (verbosity > 1) {
      std::cout << "Running Velocity Verlet integrator ..." << std::endl;
    }

    // check the validity of parameter
    checkParameters();

    totalPressure.resize(f.mesh->nVertices(), 3);
    newTotalPressure.resize(f.mesh->nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure.setZero();

    f.computeFreeEnergy();
    totalEnergy = f.E.totalE;
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
