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
  double currentNormSq;
  double pastNormSq;

  std::size_t countCG = 0;

public:
  std::size_t restartNum = 5;
  bool isBacktrack = true;
  double rho = 0.99;
  double c1 = 0.001;
  double ctol = 0.01;
  bool isAugmentedLagrangian = false;

  // std::size_t countPM = 0;

  ConjugateGradient(System &f_, double dt_, double total_time_, double tSave_,
                    double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_) {

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

// ==========================================================
// =============      FeedForward Sweep         =============
// ==========================================================
/**
 * @brief Conjugate Gradient based feedforward parameter sweep
 * @param H__, vector of sweep parameters (H0)
 * @param VP__, vector of sweep parameters (V/cam)
 * @return
 */
class DLL_PUBLIC FeedForwardSweep : public ConjugateGradient {
public:
  std::vector<double> H_;
  std::vector<double> VP_;
  FeedForwardSweep(System &f_, double dt_, double total_time_, double tSave_,
                   double tolerance_, std::string outputDir_,
                   bool isAdaptiveStep_, std::string trajFileName_,
                   std::size_t verbosity_, bool isBacktrack_, double rho_,
                   double c1_, double ctol_, bool isAugmentedLagrangian_,
                   std::size_t restartNum_, std::vector<double> H__,
                   std::vector<double> VP__)
      : ConjugateGradient(f_, dt_, total_time_, tSave_, tolerance_, outputDir_),
        H_(H__), VP_(VP__) {
    mem3dg_runtime_error(
        "FeedForwardSweep is currently not tested and maintained!");
  }
  void sweep();
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
