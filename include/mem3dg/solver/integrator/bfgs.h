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
// ==========================================================
// =============              BFGS              =============
// ==========================================================
/**
 * @brief BFGS optimizer
 * @param ctol, tolerance for termination (contraints)
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @return Success, if simulation is sucessful
 */
class DLL_PUBLIC BFGS : public Integrator {
private:
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s_protein;

public:
  const bool isBacktrack;
  const double rho;
  const double c1;
  const double ctol;
  const bool isAugmentedLagrangian;
  bool ifRestart;
  BFGS(System &system_, double characteristicTimeStep_, double totalTime_,
       double savePeriod_, double tolerance_, std::string outputDirectory_,
       bool isAdaptiveStep_, std::string trajFileName_, std::size_t verbosity_,
       bool isBacktrack_, double rho_, double c1_, double ctol_,
       bool isAugmentedLagrangian_)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_), ctol(ctol_),
        isAugmentedLagrangian(isAugmentedLagrangian_), ifRestart(false) {

    mem3dg_runtime_error("BFGS is currently not tested and maintained!");

    // print to console
    std::cout << "Running BFGS propagator ..." << std::endl;

    hess_inv.resize(system.mesh->nVertices() * 3, system.mesh->nVertices() * 3);
    hess_inv.setIdentity();
    pastPhysicalForce.resize(system.mesh->nVertices() * 3, 1);
    pastPhysicalForce.setZero();
    s.resize(system.mesh->nVertices() * 3, 1);
    s.setZero();

    hess_inv_protein.resize(system.mesh->nVertices(), system.mesh->nVertices());
    hess_inv_protein.setIdentity();
    pastPhysicalForce_protein.resize(system.mesh->nVertices(), 1);
    pastPhysicalForce_protein.setZero();
    s_protein.resize(system.mesh->nVertices(), 1);
    s_protein.setZero();

    // check the validity of parameter
    checkParameters();
  }

  /**
   * @brief BFGS function
   */
  bool integrate() override;

  /**
   * @brief BFGS stepper
   */
  void march() override;

  /**
   * @brief BFGS status computation and thresholding
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
