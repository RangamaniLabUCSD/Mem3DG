// /*
//  * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
//  *
//  * 3D Monte Carlo (MCMC) Energy Minimization Integrator.
//  *
//  * Copyright 2024- The Mem3DG Authors
//  */

// #include <Eigen/Core>
// #include <iostream>
// #include <random>

// #include <geometrycentral/surface/halfedge_mesh.h>
// #include <geometrycentral/surface/meshio.h>
// #include <geometrycentral/surface/vertex_position_geometry.h>
// #include <geometrycentral/utilities/eigen_interop_helpers.h>
// #include <geometrycentral/utilities/vector3.h>

// #include "mem3dg/meshops.h"
// #include "mem3dg/solver/integrator/mcmc.h"
// #include "mem3dg/solver/integrator/integrator.h"
// #include "mem3dg/solver/system.h"
// #include "mem3dg/type_utilities.h"

// namespace mem3dg {
// namespace solver {
// namespace integrator {
// namespace gc = ::geometrycentral;

// bool MCMC::integrate() {
//   if (ifDisableIntegrate)
//     mem3dg_runtime_error("integrate() is disabled for current construction!");

//   signal(SIGINT, signalHandler);

//   std::size_t acceptedMoves = 0;
//   std::size_t rejectedMoves = 0;

//   for (std::size_t step = 0; step < totalSteps; ++step) {
//     proposeMove();

//     // Save data periodically
//     if (step % 1000 == 0) {
//       saveData(ifOutputTrajFile, ifOutputMeshFile, ifPrintToConsole);
//       std::cout << "Step: " << step << ", Accepted: " << acceptedMoves
//                 << ", Rejected: " << rejectedMoves << std::endl;
//     }

//     // Check if termination conditions are met
//     status();
//     if (EXIT) break;
//   }

//   return SUCCESS;
// }

// void MCMC::proposeMove() {
//   static std::mt19937 rng(std::random_device{}());
//   static std::uniform_real_distribution<double> rand_displacement(-stepSize, stepSize);

//   // Store original positions
//   gcs::VertexData<gc::Vector3> originalPositions = system.geometry.vpg->inputVertexPositions;

//   // Random perturbation in 3D space
//   for (auto v : system.geometry.vpg) {
//     for (int j = 0; j < 3; ++j) {  // Move in x, y, z
//       system.geometry.vpg->inputVertexPositions(i, j) += rand_displacement(rng);
//     }
//   }

//   // Compute energy change
//   double energyOld = system.computeTotalEnergy();
//   double energyNew = system.computeTotalEnergy();
//   double deltaE = energyNew - energyOld;

//   // Accept or reject move
//   if (!metropolisAcceptance(deltaE)) {
//     system.geometry.vpg->inputVertexPositions = originalPositions;  // Revert
//   }
// }

// bool MCMC::metropolisAcceptance(double deltaE) {
//   static std::mt19937 rng(std::random_device{}());
//   static std::uniform_real_distribution<double> rand_accept(0.0, 1.0);

//   if (deltaE < 0) {
//     return true; // Always accept if energy decreases
//   }

//   double acceptanceProbability = exp(-deltaE / temperature);
//   return rand_accept(rng) < acceptanceProbability;
// }

// void MCMC::status() {
//   double energy = system.computeTotalEnergy();

//   // Exit if energy change is below tolerance
//   if (std::abs(system.mechErrorNorm) < tolerance) {
//     if (ifPrintToConsole)
//       mem3dg_print("Energy change below tolerance.");
//     EXIT = true;
//   }

//   // Exit if maximum steps reached
//   if (system.time >= totalSteps) {
//     if (ifPrintToConsole)
//       std::cout << "\nReached maximum MCMC steps." << std::endl;
//     EXIT = true;
//   }
// }

// void MCMC::checkParameters() {
//   if (temperature <= 0) {
//     mem3dg_runtime_error("Temperature must be positive for MCMC!");
//   }
//   if (stepSize <= 0) {
//     mem3dg_runtime_error("Step size must be positive for MCMC!");
//   }
// }

// } // namespace integrator
// } // namespace solver
// } // namespace mem3dg
