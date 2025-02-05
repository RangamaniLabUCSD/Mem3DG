// /*
//  * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
//  *
//  * 3D Monte Carlo (MCMC) Energy Minimization Integrator.
//  *
//  * Copyright 2024- The Mem3DG Authors
//  */

// #pragma once

// #include "mem3dg/solver/integrator/integrator.h"
// #include "mem3dg/solver/system.h"
// #include <random>

// namespace mem3dg {
// namespace solver {
// namespace integrator {

// /**
//  * @brief Monte Carlo (MCMC) Energy Minimization in 3D.
//  * @param system_ System object to be optimized.
//  * @param totalSteps_ Number of Monte Carlo iterations.
//  * @param temperature_ Metropolis temperature (controls acceptance probability).
//  * @param stepSize_ Maximum random displacement for perturbations.
//  * @param tolerance_ Energy convergence tolerance.
//  * @param outputDirectory_ Path to output directory.
//  */
// class DLL_PUBLIC MCMC : public Integrator {
// public:
//   double temperature;
//   double stepSize;
//   std::size_t totalSteps;

//   MCMC(System &system_, std::size_t totalSteps_, double temperature_, double stepSize_,
//        double tolerance_, std::string outputDirectory_)
//       : Integrator(system_, 0.0, tolerance_, outputDirectory_),
//         temperature(temperature_), stepSize(stepSize_), totalSteps(totalSteps_) {}

//   /**
//    * @brief Runs the MCMC minimization algorithm.
//    */
//   bool integrate() override;

//   /**
//    * @brief Propose a random move in 3D and decide acceptance.
//    */
//   void proposeMove();

//   /**
//    * @brief Compute acceptance probability and update the state.
//    */
//   bool metropolisAcceptance(double deltaE);

//   /**
//    * @brief Check status and convergence.
//    */
//   void status() override;

//   /**
//    * @brief Check parameters for MCMC minimization.
//    */
//   void checkParameters() override;
// };

// } // namespace integrator
// } // namespace solver
// } // namespace mem3dg
