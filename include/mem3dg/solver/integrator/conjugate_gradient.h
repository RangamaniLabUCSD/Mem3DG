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
 * @brief Conjugate Gradient propagator
 * @param ctol, tolerance for termination (contraints)
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @return Success, if simulation is sucessful
 */
class DLL_PUBLIC ConjugateGradient : public Integrator {
private:
  double currentNormSquared;
  double pastNormSquared;

  std::size_t countCG = 0;

public:
  std::size_t restartPeriod = 5;
  bool isBacktrack = true;
  double rho = 0.9;
  double c1 = 0.0005;
  double constraintTolerance = 0.01;
  bool isAugmentedLagrangian = false;

  // std::size_t countPM = 0;

  ConjugateGradient(System &system_, double characteristicTimeStep_,
                    double totalTime_, double savePeriod_, double tolerance_,
                    std::string outputDirectory_)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_) {

    // print to console
    std::cout << "Running Conjugate Gradient propagator ..." << std::endl;

    // check the validity of parameter
    checkParameters();
  }

  /**
   * @brief Conjugate Gradient driver function
   */
  bool integrate() override;

  /**
   * @brief Conjugate Gradient stepper
   */
  void march() override;

  /**
   * @brief Conjugate Gradient status computation and thresholding
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
