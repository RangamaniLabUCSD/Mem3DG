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
 * @param system_ System object to be integrated 
 * @param characteristicTimeStep_ time step, or the initial time step for backtracking algorithm 
 * @param totalTime_ total simulation time 
 * @param savePeriod_ period of saving output data (in the unit of characteristicTimeStep_) 
 * @param tolerance_ tolerance of force L2 norm for termination 
 * @param outputDirectory_ path to the output directory  
 * @param frame_ frame index, if nonzero, enable continuation mode
 */
class DLL_PUBLIC ConjugateGradient : public Integrator {
private:
  double currentNormSquared;
  double pastNormSquared;
  /// Normalized area difference to reference mesh
  double areaDifference;
  /// Normalized volume/osmotic pressure difference
  double volumeDifference;

  std::size_t countCG = 0;

public:
  std::size_t restartPeriod = 5;
  double constraintTolerance = 0.01;
  bool isAugmentedLagrangian = false;

  // std::size_t countPM = 0;

  ConjugateGradient(System &system_, double characteristicTimeStep_,
                    double totalTime_, double savePeriod_, double tolerance_,
                    std::string outputDirectory_, std::size_t frame_)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_, frame_) {
    // Initialize geometry constraints
    areaDifference = std::numeric_limits<double>::infinity();
    volumeDifference = std::numeric_limits<double>::infinity();

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

  /**
   * @brief Thresholding when adopting reduced volume parametrization
   * @param EXIT, reference to the exit flag
   * @param isAugmentedLagrangian, whether using augmented lagrangian method
   * @param dArea, normalized area difference
   * @param dVolume, normalized volume difference
   * @param ctol, exit criterion for constraint
   * @param increment, increment coefficient of penalty when using incremental
   * penalty method
   * @return
   */
  void reducedVolumeThreshold(bool &EXIT, const bool isAugmentedLagrangian,
                              const double dArea, const double dVolume,
                              const double ctol, double increment);
  /**
   * @brief Thresholding when adopting ambient pressure constraint
   * @param EXIT, reference to the exit flag
   * @param isAugmentedLagrangian, whether using augmented lagrangian method
   * @param dArea, normalized area difference
   * @param ctol, exit criterion for constraint
   * @param increment, increment coefficient of penalty when using incremental
   * penalty method
   * @return
   */
  void pressureConstraintThreshold(bool &EXIT, const bool isAugmentedLagrangian,
                                   const double dArea, const double ctol,
                                   double increment);
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
