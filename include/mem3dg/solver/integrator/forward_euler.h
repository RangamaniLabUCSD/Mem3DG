/*
 * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
 *
 * Copyright 2021- The Mem3DG Authors
 * and the project initiators Cuncheng Zhu, Christopher T. Lee, and
 * Padmini Rangamani.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Please help us support Mem3DG development by citing the research
 * papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
 * for more information.
 */

#pragma once

#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
/**
 * @brief Euler (gradient descent) time Integration
 * @param system_ System object to be integrated
 * @param characteristicTimeStep_ time step, or the initial time step for
 * backtracking algorithm
 * @param totalTime_ total simulation time
 * @param savePeriod_ period of saving output data (in the unit of
 * characteristicTimeStep_)
 * @param tolerance_ tolerance of force L2 norm for termination
 * @param outputDirectory_ path to the output directory
 * @param frame_ frame index, if nonzero, enable continuation mode
 */
class DLL_PUBLIC Euler : public Integrator {
public:
  Euler(System &system_, double characteristicTimeStep_, double totalTime_,
        double savePeriod_, double tolerance_, std::string outputDirectory_,
        std::size_t frame_ = 0)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_, frame_) {
    // check the validity of parameter
    checkParameters();
  }
  Euler(System &system_, double characteristicTimeStep_, double tolerance_,
        std::string outputDirectory_)
      : Integrator(system_, characteristicTimeStep_, tolerance_,
                   outputDirectory_) {
    checkParameters();
  }

  /**
   * @brief Forward Euler driver function
   */
  bool integrate() override;

  /**
   * @brief Forward Euler stepper
   */
  void march() override;

  // std::tuple<EigenVectorX3dr, EigenVectorX1d>
  // flowMap(EigenVectorX3dr &position, EigenVectorX1d &protein, double h);

  /**
   * @brief Forward Euler status computation and thresholding
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
