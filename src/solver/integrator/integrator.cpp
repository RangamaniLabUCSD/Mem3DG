// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2020:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/type_utilities.h"
#include "mem3dg/version.h"

#include <cmath>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace mem3dg {
namespace solver {
namespace integrator {

double Integrator::updateAdaptiveCharacteristicStep() {
  double currentMinimumSize = system.vpg->edgeLengths.raw().minCoeff();
  double currentMaximumForce =
      system.parameters.variation.isShapeVariation
          ? toMatrix(system.forces.mechanicalForce).cwiseAbs().maxCoeff()
          : toMatrix(system.forces.chemicalPotential).cwiseAbs().maxCoeff();

  double dt = (dt_size2_ratio * currentMinimumSize * currentMinimumSize) *
              (initialMaximumForce / currentMaximumForce);

  if (characteristicTimeStep / dt > 1e3) {
    mem3dg_runtime_message("Time step too small! May consider restarting the "
                           "simulation in small time scale");
    std::cout << "Current size / initial size = "
              << currentMinimumSize /
                     pow(characteristicTimeStep / dt_size2_ratio, 0.5)
              << std::endl;
    std::cout << "Current force / initial force = "
              << currentMaximumForce / initialMaximumForce << std::endl;
    EXIT = true;
    SUCCESS = false;
  }
  return dt;
}

double Integrator::backtrack(
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &&chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double positionProjection = 0;
  double chemicalProjection = 0;
  if (system.parameters.variation.isShapeVariation) {
    positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                          positionDirection.array())
                             .sum();
    if (positionProjection < 0) {
      std::cout << "\nBacktracking line search: positional velocity on uphill "
                   "direction, use bare "
                   "gradient! \n"
                << std::endl;
      positionDirection = toMatrix(system.forces.mechanicalForceVec);
      positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                            positionDirection.array())
                               .sum();
    }
  }
  if (system.parameters.variation.isProteinVariation) {
    chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                          chemicalDirection.array())
                             .sum();
    if (chemicalProjection < 0) {
      std::cout << "\nBacktracking line search: chemical direction on "
                   "uphill direction, "
                   "use bare "
                   "gradient! \n"
                << std::endl;
      chemicalDirection = system.forces.chemicalPotential.raw();
      chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                            chemicalDirection.array())
                               .sum();
    }
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  if (system.parameters.variation.isShapeVariation) {
    toMatrix(system.vpg->inputVertexPositions) += alpha * positionDirection;
  }
  if (system.parameters.variation.isProteinVariation) {
    system.proteinDensity.raw() += alpha * chemicalDirection;
  }
  system.time += alpha;
  system.updateConfigurations(false);
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (system.energy.potentialEnergy <
        (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
         c1 * alpha * (positionProjection + chemicalProjection))) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message("line search failure! Simulation "
                             "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                               toMatrix(initial_protein), previousE, true);
      std::cout << "\nError backtrace using characteristicTimeStep: \n"
                << std::endl;
      lineSearchErrorBacktrace(characteristicTimeStep,
                               toMatrix(system.vpg->inputVertexPositions),
                               toMatrix(initial_protein), previousE, true);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    if (system.parameters.variation.isShapeVariation) {
      toMatrix(system.vpg->inputVertexPositions) =
          toMatrix(initial_pos) + alpha * positionDirection;
    }
    if (system.parameters.variation.isProteinVariation) {
      system.proteinDensity.raw() =
          toMatrix(initial_protein) + alpha * chemicalDirection;
    }
    system.time = init_time + alpha;
    system.updateConfigurations(false);
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && verbosity > 3) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
    // std::cout << "mech norm: " << system.mechErrorNorm << std::endl;
    // std::cout << "chem norm: " << system.chemErrorNorm << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                             toMatrix(initial_protein), previousE, isDebug);
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations(false);
  system.computePotentialEnergy();
  return alpha;
}
double Integrator::chemicalBacktrack(
    Eigen::Matrix<double, Eigen::Dynamic, 1> &&chemicalDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double chemicalProjection = 0;
  chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                        chemicalDirection.array())
                           .sum();
  if (chemicalProjection < 0) {
    std::cout << "\nchemicalBacktracking line search: chemical direction on "
                 "uphill direction, "
                 "use bare "
                 "gradient! \n"
              << std::endl;
    chemicalDirection = system.forces.chemicalPotential.raw();
    chemicalProjection = (system.forces.chemicalPotential.raw().array() *
                          chemicalDirection.array())
                             .sum();
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  system.proteinDensity.raw() += alpha * chemicalDirection;
  system.time += alpha;
  system.updateConfigurations(false);
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (system.energy.potentialEnergy <
        (previousE.potentialEnergy - c1 * alpha * chemicalProjection)) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message(
          "\nchemicalBacktrack: line search failure! Simulation "
          "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                               toMatrix(initial_protein), previousE, true);
      std::cout << "\nError backtrace using characteristicTimeStep: \n"
                << std::endl;
      lineSearchErrorBacktrace(characteristicTimeStep, toMatrix(initial_pos),
                               toMatrix(initial_protein), previousE, true);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    system.proteinDensity.raw() =
        toMatrix(initial_protein) + alpha * chemicalDirection;
    system.time = init_time + alpha;
    system.updateConfigurations(false);
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && verbosity > 3) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    std::cout << "\nchemicalBacktrack: debugging \n" << std::endl;
    lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                             toMatrix(initial_protein), previousE, isDebug);
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations(false);
  system.computePotentialEnergy();
  return alpha;
}

double Integrator::mechanicalBacktrack(
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection, double rho,
    double c1) {

  // cache energy of the last time step
  const Energy previousE = system.energy;

  // validate the directions
  double positionProjection = 0;
  positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                        positionDirection.array())
                           .sum();
  if (positionProjection < 0) {
    std::cout
        << "\nmechanicalBacktrack line search: positional velocity on uphill "
           "direction, use bare "
           "gradient! \n"
        << std::endl;
    positionDirection = toMatrix(system.forces.mechanicalForceVec);
    positionProjection = (toMatrix(system.forces.mechanicalForceVec).array() *
                          positionDirection.array())
                             .sum();
  }

  // calculate initial energy as reference level
  gc::VertexData<gc::Vector3> initial_pos(*system.mesh);
  initial_pos = system.vpg->inputVertexPositions;
  gc::VertexData<double> initial_protein(*system.mesh,
                                         system.proteinDensity.raw());
  const double init_time = system.time;

  // declare variables used in backtracking iterations
  double alpha = characteristicTimeStep;
  std::size_t count = 0;

  // zeroth iteration
  toMatrix(system.vpg->inputVertexPositions) += alpha * positionDirection;
  system.time += alpha;
  system.updateConfigurations(false);
  system.computePotentialEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if ((system.energy.potentialEnergy <
         (previousE.potentialEnergy + system.computeIntegratedPower(alpha) -
          c1 * alpha * positionProjection)) &&
        std::isfinite(system.energy.potentialEnergy)) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * characteristicTimeStep) {
      mem3dg_runtime_message(
          "\nmechanicalBacktrack: line search failure! Simulation "
          "stopped. \n");
      std::cout << "\nError backtrace using alpha: \n" << std::endl;
      lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                               toMatrix(initial_protein), previousE, true);
      std::cout << "\nError backtrace using characterisiticTimeStep: \n"
                << std::endl;
      lineSearchErrorBacktrace(characteristicTimeStep, toMatrix(initial_pos),
                               toMatrix(initial_protein), previousE, true);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    toMatrix(system.vpg->inputVertexPositions) =
        toMatrix(initial_pos) + alpha * positionDirection;

    system.time = init_time + alpha;
    system.updateConfigurations(false);
    system.computePotentialEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != characteristicTimeStep && verbosity > 3) {
    std::cout << "alpha: " << characteristicTimeStep << " -> " << alpha
              << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    std::cout << "\nmechanicalBacktrack: debugging \n" << std::endl;
    lineSearchErrorBacktrace(alpha, toMatrix(initial_pos),
                             toMatrix(initial_protein), previousE, isDebug);
  }

  // recover the initial configuration
  system.time = init_time;
  system.proteinDensity = initial_protein;
  system.vpg->inputVertexPositions = initial_pos;
  system.updateConfigurations(false);
  system.computePotentialEnergy();
  return alpha;
}

void Integrator::lineSearchErrorBacktrace(
    const double alpha, const EigenVectorX3dr currentPosition,
    const EigenVectorX1d currentProteinDensity, const Energy previousEnergy,
    bool runAll) {
  std::cout << "\nlineSearchErrorBacktracking ..." << std::endl;

  // cache the energy when applied the total force
  // system.proteinDensity.raw() = currentProteinDensity;
  // toMatrix(system.vpg->inputVertexPositions) =
  //     currentPosition +
  //     alpha * system.forces.maskForce(
  //                 toMatrix(system.forces.osmoticForceVec) +
  //                 toMatrix(system.forces.capillaryForceVec) +
  //                 toMatrix(system.forces.bendingForceVec) +
  //                 toMatrix(system.forces.lineCapillaryForceVec) +
  //                 toMatrix(system.forces.adsorptionForceVec) +
  //                 toMatrix(system.forces.aggregationForceVec) +
  //                 toMatrix(system.forces.externalForceVec) +
  //                 toMatrix(system.forces.selfAvoidanceForceVec));
  // // * toMatrix(system.forces.mechanicalForceVec);
  // system.updateConfigurations(false);
  if (system.parameters.external.Kf != 0)
    system.computeExternalWork(system.time, alpha);
  system.computeTotalEnergy();
  const Energy totalForceEnergy{system.energy};

  // test if total potential energy increases
  if (runAll ||
      totalForceEnergy.potentialEnergy > previousEnergy.potentialEnergy) {

    // report the finding
    std::cout << "\nWith F_tol, potE has increased "
              << totalForceEnergy.potentialEnergy -
                     previousEnergy.potentialEnergy
              << " from " << previousEnergy.potentialEnergy << " to "
              << totalForceEnergy.potentialEnergy << std::endl;

    // test if bending energy increases
    if (runAll ||
        totalForceEnergy.bendingEnergy > previousEnergy.bendingEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, BE has increased "
                << totalForceEnergy.bendingEnergy - previousEnergy.bendingEnergy
                << " from " << previousEnergy.bendingEnergy << " to "
                << totalForceEnergy.bendingEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.bendingForceVec));
      system.updateConfigurations(false);

      // test if bending energy increases
      system.computeBendingEnergy();
      if (runAll ||
          system.energy.bendingEnergy > previousEnergy.bendingEnergy) {
        std::cout << "With only bending force, BE has increased "
                  << system.energy.bendingEnergy - previousEnergy.bendingEnergy
                  << " from " << previousEnergy.bendingEnergy << " to "
                  << system.energy.bendingEnergy << ", expected dBE: "
                  << -alpha *
                         system.forces
                             .maskForce(toMatrix(system.forces.bendingForceVec))
                             .squaredNorm()
                  << std::endl;
      }

      // perturb the configuration
      toMatrix(system.vpg->inputVertexPositions) = currentPosition;
      system.proteinDensity.raw() =
          currentProteinDensity +
          alpha * system.parameters.proteinMobility *
              system.forces.maskProtein(system.forces.bendingPotential.raw());
      system.updateConfigurations(false);

      // test if bending energy increases
      system.computeBendingEnergy();
      if (runAll ||
          system.energy.bendingEnergy > previousEnergy.bendingEnergy) {
        std::cout << "With only bending potential, BE has increased "
                  << system.energy.bendingEnergy - previousEnergy.bendingEnergy
                  << " from " << previousEnergy.bendingEnergy << " to "
                  << system.energy.bendingEnergy << ", expected dBE: "
                  << -alpha * system.parameters.proteinMobility *
                         system.forces
                             .maskProtein(system.forces.bendingPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }

    // test if deviatoric energy increases
    if (runAll ||
        totalForceEnergy.deviatoricEnergy > previousEnergy.deviatoricEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, DevE has increased "
                << totalForceEnergy.deviatoricEnergy -
                       previousEnergy.deviatoricEnergy
                << " from " << previousEnergy.deviatoricEnergy << " to "
                << totalForceEnergy.deviatoricEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.deviatoricForceVec));
      system.updateConfigurations(false);

      // test if deviatoric energy increases
      system.computeDeviatoricEnergy();
      if (runAll ||
          system.energy.deviatoricEnergy > previousEnergy.deviatoricEnergy) {
        std::cout << "With only deviatoric force, DevE has increased "
                  << system.energy.deviatoricEnergy -
                         previousEnergy.deviatoricEnergy
                  << " from " << previousEnergy.deviatoricEnergy << " to "
                  << system.energy.deviatoricEnergy << ", expected dDevE: "
                  << -alpha * system.forces
                                  .maskForce(toMatrix(
                                      system.forces.deviatoricForceVec))
                                  .squaredNorm()
                  << std::endl;
      }

      // perturb the configuration
      toMatrix(system.vpg->inputVertexPositions) = currentPosition;
      system.proteinDensity.raw() =
          currentProteinDensity +
          alpha * system.parameters.proteinMobility *
              system.forces.maskProtein(
                  system.forces.deviatoricPotential.raw());
      system.updateConfigurations(false);

      // test if deviatoric energy increases
      system.computeDeviatoricEnergy();
      if (runAll ||
          system.energy.deviatoricEnergy > previousEnergy.deviatoricEnergy) {
        std::cout << "With only deviatoric potential, DevE has increased "
                  << system.energy.deviatoricEnergy -
                         previousEnergy.deviatoricEnergy
                  << " from " << previousEnergy.deviatoricEnergy << " to "
                  << system.energy.deviatoricEnergy << ", expected dDevE: "
                  << -alpha * system.parameters.proteinMobility *
                         system.forces
                             .maskProtein(
                                 system.forces.deviatoricPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }

    // test if surface energy increases
    if (runAll ||
        totalForceEnergy.surfaceEnergy > previousEnergy.surfaceEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, sE has increased "
                << totalForceEnergy.surfaceEnergy - previousEnergy.surfaceEnergy
                << " from " << previousEnergy.surfaceEnergy << " to "
                << totalForceEnergy.surfaceEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.capillaryForceVec));
      system.updateConfigurations(false);
      system.computeSurfaceEnergy();
      if (runAll ||
          system.energy.surfaceEnergy > previousEnergy.surfaceEnergy) {
        std::cout << "With only capillary force, sE has increased "
                  << system.energy.surfaceEnergy - previousEnergy.surfaceEnergy
                  << " from " << previousEnergy.surfaceEnergy << " to "
                  << system.energy.surfaceEnergy << ", expected dsE: "
                  << -alpha * system.forces
                                  .maskForce(
                                      toMatrix(system.forces.capillaryForceVec))
                                  .squaredNorm()
                  << std::endl;
      }
    }

    // test if pressure energy increases
    if (runAll ||
        totalForceEnergy.pressureEnergy > previousEnergy.pressureEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, pE has increased "
                << totalForceEnergy.pressureEnergy -
                       previousEnergy.pressureEnergy
                << " from " << previousEnergy.pressureEnergy << " to "
                << totalForceEnergy.pressureEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.osmoticForceVec));
      system.updateConfigurations(false);
      system.computePressureEnergy();
      if (runAll ||
          system.energy.pressureEnergy > previousEnergy.pressureEnergy) {
        std::cout << "With only osmotic force, pE has increased "
                  << system.energy.pressureEnergy -
                         previousEnergy.pressureEnergy
                  << " from " << previousEnergy.pressureEnergy << " to "
                  << system.energy.pressureEnergy << ", expected dpE: "
                  << -alpha *
                         system.forces
                             .maskForce(toMatrix(system.forces.osmoticForceVec))
                             .squaredNorm()
                  << std::endl;
      }
    }

    // test if adsorption energy increases
    if (runAll ||
        totalForceEnergy.adsorptionEnergy > previousEnergy.adsorptionEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, aE has increased "
                << totalForceEnergy.adsorptionEnergy -
                       previousEnergy.adsorptionEnergy
                << " from " << previousEnergy.adsorptionEnergy << " to "
                << totalForceEnergy.adsorptionEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.adsorptionForceVec));
      system.updateConfigurations(false);
      system.computeAdsorptionEnergy();
      if (runAll ||
          system.energy.adsorptionEnergy > previousEnergy.adsorptionEnergy) {
        std::cout << "With only adsorption force, aE has increased "
                  << system.energy.adsorptionEnergy -
                         previousEnergy.adsorptionEnergy
                  << " from " << previousEnergy.adsorptionEnergy << " to "
                  << system.energy.adsorptionEnergy << ", expected daE: "
                  << -alpha * system.forces
                                  .maskForce(toMatrix(
                                      system.forces.adsorptionForceVec))
                                  .squaredNorm()
                  << std::endl;
      }

      // test single-force-energy computation
      // perturb the configuration
      toMatrix(system.vpg->inputVertexPositions) = currentPosition;
      system.proteinDensity.raw() =
          currentProteinDensity +
          alpha * system.parameters.proteinMobility *
              system.forces.maskProtein(
                  system.forces.adsorptionPotential.raw());
      system.updateConfigurations(false);
      system.computeAdsorptionEnergy();
      if (runAll ||
          system.energy.adsorptionEnergy > previousEnergy.adsorptionEnergy) {
        std::cout << "With only adsorption potential, aE has increased "
                  << system.energy.adsorptionEnergy -
                         previousEnergy.adsorptionEnergy
                  << " from " << previousEnergy.adsorptionEnergy << " to "
                  << system.energy.adsorptionEnergy << ", expected daE: "
                  << -alpha * system.parameters.proteinMobility *
                         system.forces
                             .maskProtein(
                                 system.forces.adsorptionPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }

    // test if aggregation energy increases
    if (runAll ||
        totalForceEnergy.aggregationEnergy > previousEnergy.aggregationEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, aggE has increased "
                << totalForceEnergy.aggregationEnergy -
                       previousEnergy.aggregationEnergy
                << " from " << previousEnergy.aggregationEnergy << " to "
                << totalForceEnergy.aggregationEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.aggregationForceVec));
      system.updateConfigurations(false);
      system.computeAggregationEnergy();
      if (runAll ||
          system.energy.aggregationEnergy > previousEnergy.aggregationEnergy) {
        std::cout << "With only aggregation force, aggE has increased "
                  << system.energy.aggregationEnergy -
                         previousEnergy.aggregationEnergy
                  << " from " << previousEnergy.aggregationEnergy << " to "
                  << system.energy.aggregationEnergy << ", expected daggE: "
                  << -alpha * system.forces
                                  .maskForce(toMatrix(
                                      system.forces.aggregationForceVec))
                                  .squaredNorm()
                  << std::endl;
      }

      // test single-force-energy computation
      // perturb the configuration
      toMatrix(system.vpg->inputVertexPositions) = currentPosition;
      system.proteinDensity.raw() =
          currentProteinDensity +
          alpha * system.parameters.proteinMobility *
              system.forces.maskProtein(
                  system.forces.aggregationPotential.raw());
      system.updateConfigurations(false);
      system.computeAggregationEnergy();
      if (runAll ||
          system.energy.aggregationEnergy > previousEnergy.aggregationEnergy) {
        std::cout << "With only aggregation potential, aggE has increased "
                  << system.energy.aggregationEnergy -
                         previousEnergy.aggregationEnergy
                  << " from " << previousEnergy.aggregationEnergy << " to "
                  << system.energy.aggregationEnergy << ", expected daggE: "
                  << -alpha * system.parameters.proteinMobility *
                         system.forces
                             .maskProtein(
                                 system.forces.aggregationPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }

    // test if dirichlet energy increases
    if (runAll ||
        system.energy.dirichletEnergy > previousEnergy.dirichletEnergy) {

      // report the finding
      std::cout << "\nWith F_tol, dE has increased "
                << system.energy.dirichletEnergy -
                       previousEnergy.dirichletEnergy
                << " from " << previousEnergy.dirichletEnergy << " to "
                << system.energy.dirichletEnergy << std::endl;

      // test single-force-energy computation
      // perturb the configuration
      system.proteinDensity.raw() = currentProteinDensity;
      toMatrix(system.vpg->inputVertexPositions) =
          currentPosition + alpha * system.forces.maskForce(toMatrix(
                                        system.forces.lineCapillaryForceVec));
      system.updateConfigurations(false);
      system.computeDirichletEnergy();
      if (runAll ||
          system.energy.dirichletEnergy > previousEnergy.dirichletEnergy) {
        std::cout << "With only line tension force, dE has increased "
                  << system.energy.dirichletEnergy -
                         previousEnergy.dirichletEnergy
                  << " from " << previousEnergy.dirichletEnergy << " to "
                  << system.energy.dirichletEnergy << ", expected ddE: "
                  << -alpha * system.forces
                                  .maskForce(toMatrix(
                                      system.forces.lineCapillaryForceVec))
                                  .squaredNorm()
                  << std::endl;
      }

      // test single-force-energy computation
      // perturb the configuration
      toMatrix(system.vpg->inputVertexPositions) = currentPosition;
      system.proteinDensity.raw() =
          currentProteinDensity +
          alpha * system.parameters.proteinMobility *
              system.forces.maskProtein(system.forces.diffusionPotential.raw());
      system.updateConfigurations(false);
      system.computeDirichletEnergy();
      if (runAll ||
          system.energy.dirichletEnergy > previousEnergy.dirichletEnergy) {
        std::cout << "With only diffusion potential, dE has increased "
                  << system.energy.dirichletEnergy -
                         previousEnergy.dirichletEnergy
                  << " from " << previousEnergy.dirichletEnergy << " to "
                  << system.energy.dirichletEnergy << ", expected ddE: "
                  << -alpha * system.parameters.proteinMobility *
                         system.forces
                             .maskProtein(
                                 system.forces.diffusionPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }
  }

  // test if self avoidance energy increases
  if (runAll || system.energy.selfAvoidancePenalty >
                    previousEnergy.selfAvoidancePenalty) {

    // report the finding
    std::cout << "\nWith F_tol, selfE has increased "
              << system.energy.selfAvoidancePenalty -
                     previousEnergy.selfAvoidancePenalty
              << " from " << previousEnergy.selfAvoidancePenalty << " to "
              << system.energy.selfAvoidancePenalty << std::endl;

    // test single-force-energy computation
    // perturb the configuration
    system.proteinDensity.raw() = currentProteinDensity;
    toMatrix(system.vpg->inputVertexPositions) =
        currentPosition + alpha * system.forces.maskForce(toMatrix(
                                      system.forces.selfAvoidanceForceVec));
    system.updateConfigurations(false);
    system.computeSelfAvoidanceEnergy();
    if (runAll || system.energy.selfAvoidancePenalty >
                      previousEnergy.selfAvoidancePenalty) {
      std::cout
          << "With only self avoidance penalty force, selfE has increased "
          << system.energy.selfAvoidancePenalty -
                 previousEnergy.selfAvoidancePenalty
          << " from " << previousEnergy.selfAvoidancePenalty << " to "
          << system.energy.selfAvoidancePenalty << ", expected dselfE: "
          << -alpha *
                 system.forces
                     .maskForce(toMatrix(system.forces.selfAvoidanceForceVec))
                     .squaredNorm()
          << std::endl;
    }
  }

  // test if interior penalty energy increases
  if (runAll || totalForceEnergy.proteinInteriorPenalty >
                    previousEnergy.proteinInteriorPenalty) {

    // report the finding
    std::cout << "\nWith F_tol, inPE has increased "
              << totalForceEnergy.proteinInteriorPenalty -
                     previousEnergy.proteinInteriorPenalty
              << " from " << previousEnergy.proteinInteriorPenalty << " to "
              << totalForceEnergy.proteinInteriorPenalty << std::endl;

    // test single-force-energy computation
    // perturb the configuration
    toMatrix(system.vpg->inputVertexPositions) = currentPosition;
    system.proteinDensity.raw() =
        currentProteinDensity +
        alpha * system.parameters.proteinMobility *
            system.forces.maskProtein(
                system.forces.interiorPenaltyPotential.raw());
    system.updateConfigurations(false);
    system.computeProteinInteriorPenalty();
    if (runAll || system.energy.proteinInteriorPenalty >
                      previousEnergy.proteinInteriorPenalty) {
      std::cout
          << "With only protein interior penalty potential, inPE has increased "
          << system.energy.proteinInteriorPenalty -
                 previousEnergy.proteinInteriorPenalty
          << " from " << previousEnergy.proteinInteriorPenalty << " to "
          << system.energy.proteinInteriorPenalty << ", expected dinPE: "
          << -alpha * system.parameters.proteinMobility *
                 system.forces
                     .maskProtein(system.forces.interiorPenaltyPotential.raw())
                     .squaredNorm()
          << std::endl;
    }
  }

  // test if total force is doing negative work against external force field
  if (runAll || totalForceEnergy.externalWork < previousEnergy.externalWork) {
    std::cout
        << "\nF_tol is doing negative work against external force field by "
        << previousEnergy.externalWork - totalForceEnergy.externalWork
        << std::endl;
  }

  // test if total kinetic energy increases
  if (runAll || totalForceEnergy.kineticEnergy > previousEnergy.kineticEnergy) {
    std::cout << "\nWith F_tol, kE has increased "
              << totalForceEnergy.kineticEnergy - previousEnergy.kineticEnergy
              << " from " << previousEnergy.kineticEnergy << " to "
              << totalForceEnergy.kineticEnergy << std::endl;
  }
}

void Integrator::finitenessErrorBacktrace() {

  if (!std::isfinite(timeStep)) {
    EXIT = true;
    SUCCESS = false;
    mem3dg_runtime_message("time step is not finite!");
  }

  if (!std::isfinite(system.mechErrorNorm)) {
    EXIT = true;
    SUCCESS = false;

    if (!std::isfinite(toMatrix(system.velocity).norm())) {
      mem3dg_runtime_message("Velocity is not finite!");
    }

    if (!std::isfinite(toMatrix(system.forces.mechanicalForceVec).norm())) {
      if (!std::isfinite(toMatrix(system.forces.capillaryForceVec).norm())) {
        mem3dg_runtime_message("Capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.adsorptionForceVec).norm())) {
        mem3dg_runtime_message("Adsorption force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.aggregationForceVec).norm())) {
        mem3dg_runtime_message("Aggregation force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.bendingForceVec).norm())) {
        mem3dg_runtime_message("Bending force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.deviatoricForceVec).norm())) {
        mem3dg_runtime_message("Deviatoric force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.osmoticForceVec).norm())) {
        mem3dg_runtime_message("Osmotic force is not finite!");
      }
      if (!std::isfinite(
              toMatrix(system.forces.lineCapillaryForceVec).norm())) {
        mem3dg_runtime_message("Line capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.externalForceVec).norm())) {
        mem3dg_runtime_message("External force is not finite!");
      }
      if (!std::isfinite(
              toMatrix(system.forces.selfAvoidanceForceVec).norm())) {
        mem3dg_runtime_message("Self avoidance force is not finite!");
      }
    }
  }

  if (!std::isfinite(system.chemErrorNorm)) {
    EXIT = true;
    SUCCESS = false;

    if (!std::isfinite(toMatrix(system.proteinVelocity).norm())) {
      mem3dg_runtime_message("Protein velocity is not finite!");
    }

    if (!std::isfinite(toMatrix(system.forces.chemicalPotential).norm())) {
      if (!std::isfinite(toMatrix(system.forces.bendingPotential).norm())) {
        mem3dg_runtime_message("Bending Potential is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.deviatoricPotential).norm())) {
        mem3dg_runtime_message("Deviatoric Potential is not finite!");
      }
      if (!std::isfinite(
              toMatrix(system.forces.interiorPenaltyPotential).norm())) {
        mem3dg_runtime_message(
            "Protein interior penalty potential is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.diffusionPotential).norm())) {
        mem3dg_runtime_message("Diffusion potential is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.adsorptionPotential).norm())) {
        mem3dg_runtime_message("Adsorption potential is not finite!");
      }
      if (!std::isfinite(toMatrix(system.forces.aggregationPotential).norm())) {
        mem3dg_runtime_message("Aggregation potential is not finite!");
      }
    }
  }

  if (!std::isfinite(system.energy.totalEnergy)) {
    EXIT = true;
    SUCCESS = false;
    if (!std::isfinite(system.energy.kineticEnergy)) {
      mem3dg_runtime_message("Kinetic energy is not finite!");
    }
    if (!std::isfinite(system.energy.externalWork)) {
      mem3dg_runtime_message("External work is not finite!");
    }
    if (!std::isfinite(system.energy.potentialEnergy)) {
      if (!std::isfinite(system.energy.bendingEnergy)) {
        mem3dg_runtime_message("Bending energy is not finite!");
      }
      if (!std::isfinite(system.energy.deviatoricEnergy)) {
        mem3dg_runtime_message("Deviatoric energy is not finite!");
      }
      if (!std::isfinite(system.energy.surfaceEnergy)) {
        mem3dg_runtime_message("Surface energy is not finite!");
      }
      if (!std::isfinite(system.energy.pressureEnergy)) {
        mem3dg_runtime_message("Pressure energy is not finite!");
      }
      if (!std::isfinite(system.energy.adsorptionEnergy)) {
        mem3dg_runtime_message("Adsorption energy is not finite!");
      }
      if (!std::isfinite(system.energy.aggregationEnergy)) {
        mem3dg_runtime_message("Aggregation energy is not finite!");
      }
      if (!std::isfinite(system.energy.dirichletEnergy)) {
        mem3dg_runtime_message("Line tension energy is not finite!");
      }
      if (!std::isfinite(system.energy.proteinInteriorPenalty)) {
        mem3dg_runtime_message(
            "Protein interior penalty energy is not finite!");
      }
      if (!std::isfinite(system.energy.selfAvoidancePenalty)) {
        mem3dg_runtime_message(
            "Membrane self-avoidance penalty energy is not finite!");
      }
    }
  }
}

void Integrator::pressureConstraintThreshold(bool &EXIT,
                                             const bool isAugmentedLagrangian,
                                             const double dArea,
                                             const double ctol,
                                             double increment) {
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    if (isAugmentedLagrangian) { // augmented Lagrangian method
      if (dArea < ctol) {        // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG] = [" << system.parameters.tension.lambdaSG
                  << ", "
                  << "]";
        system.parameters.tension.lambdaSG +=
            system.parameters.tension.Ksg *
            (system.surfaceArea - system.parameters.tension.At) /
            system.parameters.tension.At;
        std::cout << " -> [" << system.parameters.tension.lambdaSG << "]"
                  << std::endl;
      }
    } else {              // incremental harmonic penalty method
      if (dArea < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[Ksg] = [" << system.parameters.tension.Ksg << "]";
        system.parameters.tension.Ksg *= increment;
        std::cout << " -> [" << system.parameters.tension.Ksg << "]"
                  << std::endl;
      }
    }
  }
}

void Integrator::reducedVolumeThreshold(bool &EXIT,
                                        const bool isAugmentedLagrangian,
                                        const double dArea,
                                        const double dVolume, const double ctol,
                                        double increment) {
  if (system.mechErrorNorm < tolerance && system.chemErrorNorm < tolerance) {
    if (isAugmentedLagrangian) {            // augmented Lagrangian method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n.tension[lambdaSG, lambdaV] = ["
                  << system.parameters.tension.lambdaSG << ", "
                  << system.parameters.osmotic.lambdaV << "]";
        system.parameters.tension.lambdaSG +=
            system.parameters.tension.Ksg *
            (system.surfaceArea - system.parameters.tension.At) /
            system.parameters.tension.At;
        system.parameters.osmotic.lambdaV +=
            system.parameters.osmotic.Kv *
            (system.volume - system.parameters.osmotic.Vt) /
            system.parameters.osmotic.Vt;
        std::cout << " -> [" << system.parameters.tension.lambdaSG << ", "
                  << system.parameters.osmotic.lambdaV << "]" << std::endl;
      }
    } else { // incremental harmonic penalty method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      }

      // iterate if not
      if (dArea > ctol) {
        std::cout << "\n[Ksg] = [" << system.parameters.tension.Ksg << "]";
        system.parameters.tension.Ksg *= 1.3;
        std::cout << " -> [" << system.parameters.tension.Ksg << "]"
                  << std::endl;
      }
      if (dVolume > ctol) {
        std::cout << "\n[Kv] = [" << system.parameters.osmotic.Kv << "]";
        system.parameters.osmotic.Kv *= 1.3;
        std::cout << " -> [" << system.parameters.osmotic.Kv << "]"
                  << std::endl;
      }
    }
  }
}

void Integrator::saveData() {
  // threshold of verbosity level to output ply file
  int outputPly = 0;

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (verbosity > 0) {
    // saveNetcdfData();
    saveMutableNetcdfData();
  }
  outputPly = 3;
#endif

  // save variable to richData and save ply file
  if (verbosity > outputPly) {
    char buffer[50];
    sprintf(buffer, isJustGeometryPly ? "/frame%d.obj" : "/frame%d.ply",
            (int)frame);
    system.saveRichData(outputDirectory + "/" + std::string(buffer),
                        isJustGeometryPly);
  }

  // print in-progress information in the console
  if (verbosity > 1) {
    std::cout << "\n"
              << "t: " << system.time << ", "
              << "n: " << frame << ", "
              << "isSmooth: " << system.isSmooth << "\n"
              << "dA/Area: " << areaDifference << "/" << system.surfaceArea
              << ", "
              << "dVP/Volume: " << volumeDifference << "/" << system.volume
              << ", "
              << "h: "
              << toMatrix(system.vpg->inputVertexPositions).col(2).maxCoeff()
              << "\n"
              << "E_total: " << system.energy.totalEnergy << "\n"
              << "E_kin: " << system.energy.kineticEnergy << "\n"
              << "E_pot: " << system.energy.potentialEnergy << "\n"
              << "W_ext: " << system.energy.externalWork << "\n"
              << "|e|Mech: " << system.mechErrorNorm << "\n"
              << "|e|Chem: " << system.chemErrorNorm << "\n"
              << "H: ["
              << (system.vpg->vertexMeanCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (system.vpg->vertexMeanCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]"
              << "\n"
              << "K: ["
              << (system.vpg->vertexGaussianCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (system.vpg->vertexGaussianCurvatures.raw().array() /
                  system.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]"
              << "\n"
              << "phi: [" << system.proteinDensity.raw().minCoeff() << ","
              << system.proteinDensity.raw().maxCoeff() << "]" << std::endl;
    // << "COM: "
    // << gc::EigenMap<double,
    // 3>(f.vpg->inputVertexPositions).colwise().sum() /
    //         f.vpg->inputVertexPositions.raw().rows()
    // << "\n"
  }
  // break loop if EXIT flag is on
  if (EXIT) {
    if (verbosity > 0) {
      std::cout << "Simulation " << (SUCCESS ? "finished" : "failed")
                << ", and data saved to " + outputDirectory << std::endl;
      if (verbosity > 2) {
        system.saveRichData(outputDirectory + "/out.ply");
      }
    }
  }

  frame++;
}

void Integrator::markFileName(std::string marker_str) {
  std::string dirPath = outputDirectory;

  const char *marker = marker_str.c_str();

  char *file = new char[trajFileName.size() + 1];
  std::copy(trajFileName.begin(), trajFileName.end(), file);
  file[trajFileName.size()] = '\0';

  char fileMarked[50]{"/"}, oldNC[150]{"/"}, newNC[150]{"/"};

  // sprintf(fileMarked, "/traj_H_%d_VP_%d_failed.nc", int(H * 100),
  //         int(VP * 100));

  // split the extension and file name
  const char *ext = strchr(file, '.');

  // name fileMarked to be the file name
  std::strncpy(fileMarked, file, ext - file);

  // name fileMarked to be file name + the marker + extension
  std::strcat(fileMarked, marker);
  std::strcat(fileMarked, ext);
  fileMarked[ext - file + sizeof(marker) + sizeof(ext)] = '\0';

  // append the directory path and copy to oldNC and newNC
  std::strcpy(oldNC, dirPath.c_str());
  std::strcpy(newNC, dirPath.c_str());
  std::strcat(oldNC, file);
  std::strcat(newNC, fileMarked);

  // rename file
  rename(oldNC, newNC);
  delete[] file;
}

#ifdef MEM3DG_WITH_NETCDF
void Integrator::createNetcdfFile() {
  // initialize netcdf traj file
  trajFile.createNewFile(outputDirectory + "/" + trajFileName, *system.mesh,
                         *system.vpg, TrajFile::NcFile::replace);
  trajFile.writeMask(toMatrix(system.forces.forceMask).rowwise().sum());
  if (!system.mesh->hasBoundary()) {
    trajFile.writeRefSurfArea(system.parameters.tension.At);
  }
}

void Integrator::createMutableNetcdfFile() {
  // initialize netcdf traj file
  mutableTrajFile.createNewFile(outputDirectory + "/" + trajFileName,
                                TrajFile::NcFile::replace);
  // mutableTrajFile.writeMask(toMatrix(f.forces.forceMask).rowwise().sum());
  // if (!f.mesh->hasBoundary()) {
  //   mutableTrajFile.writeRefSurfArea(f.parameters.tension.At);
  // }
}

void Integrator::saveNetcdfData() {
  std::size_t idx = trajFile.nFrames();

  // scalar quantities
  // write time
  trajFile.writeTime(idx, system.time);
  trajFile.writeIsSmooth(idx, system.isSmooth);
  // write geometry
  trajFile.writeVolume(idx, system.volume);
  trajFile.writeSurfArea(idx,
                         system.mesh->hasBoundary()
                             ? system.surfaceArea - system.parameters.tension.At
                             : system.surfaceArea);
  trajFile.writeHeight(
      idx, toMatrix(system.vpg->inputVertexPositions).col(2).maxCoeff());
  // write energies
  trajFile.writeBendEnergy(idx, system.energy.bendingEnergy);
  trajFile.writeSurfEnergy(idx, system.energy.surfaceEnergy);
  trajFile.writePressEnergy(idx, system.energy.pressureEnergy);
  trajFile.writeKineEnergy(idx, system.energy.kineticEnergy);
  trajFile.writeAdspEnergy(idx, system.energy.adsorptionEnergy);
  trajFile.writeLineEnergy(idx, system.energy.dirichletEnergy);
  trajFile.writeTotalEnergy(idx, system.energy.totalEnergy);
  // write Norms
  trajFile.writeErrorNorm(idx, system.mechErrorNorm);
  trajFile.writeChemErrorNorm(idx, system.chemErrorNorm);
  trajFile.writeBendNorm(idx,
                         system.computeNorm(system.forces.bendingForce.raw()));
  trajFile.writeSurfNorm(
      idx, system.computeNorm(system.forces.capillaryForce.raw()));
  trajFile.writePressNorm(idx,
                          system.computeNorm(system.forces.osmoticForce.raw()));
  trajFile.writeLineNorm(
      idx, system.computeNorm(system.forces.lineCapillaryForce.raw()));

  // vector quantities
  if (!system.meshProcessor.meshMutator.isSplitEdge &&
      !system.meshProcessor.meshMutator.isCollapseEdge) {
    // write velocity
    trajFile.writeVelocity(idx, toMatrix(system.velocity));
    // write protein density distribution
    trajFile.writeProteinDensity(idx, system.proteinDensity.raw());

    // write geometry
    trajFile.writeCoords(idx, toMatrix(system.vpg->inputVertexPositions));
    trajFile.writeTopoFrame(idx,
                            system.mesh->getFaceVertexMatrix<std::uint32_t>());
    trajFile.writeMeanCurvature(idx,
                                system.vpg->vertexMeanCurvatures.raw().array() /
                                    system.vpg->vertexDualAreas.raw().array());
    trajFile.writeGaussCurvature(
        idx, system.vpg->vertexGaussianCurvatures.raw().array() /
                 system.vpg->vertexDualAreas.raw().array());
    trajFile.writeSponCurvature(idx, system.H0.raw());
    // fd.writeAngles(idx, f.vpg.cornerAngles.raw());
    // fd.writeH_H0_diff(idx,
    //                   ((f.H - f.H0).array() * (f.H -
    //                   f.H0).array()).matrix());

    // write pressures
    trajFile.writeBendingForce(idx, system.forces.bendingForce.raw());
    trajFile.writeCapillaryForce(idx, system.forces.capillaryForce.raw());
    trajFile.writeLineForce(idx, system.forces.lineCapillaryForce.raw());
    trajFile.writeOsmoticForce(idx, system.forces.osmoticForce.raw());
    trajFile.writeExternalForce(idx, system.forces.externalForce.raw());
    trajFile.writePhysicalForce(idx, system.forces.mechanicalForce.raw());
    trajFile.writeChemicalPotential(idx, system.forces.chemicalPotential.raw());
  }
}

void Integrator::saveMutableNetcdfData() {
  std::size_t idx = mutableTrajFile.nFrames();

  // scalar quantities
  // write time
  mutableTrajFile.writeTime(idx, system.time);

  // write dynamic properties
  mutableTrajFile.writeVelocity(idx, system.velocity);
  if (system.parameters.external.Kf != 0)
    mutableTrajFile.writeExternalForce(idx, system.forces.externalForceVec);

  // write static properties
  mutableTrajFile.writeCoords(idx, *system.vpg);
  mutableTrajFile.writeTopology(idx, *system.mesh);
  mutableTrajFile.writeProteinDensity(idx, system.proteinDensity);
  mutableTrajFile.sync();
}
#endif

void Integrator::getParameterLog(std::string inputMesh) {
  std::ofstream myfile(outputDirectory + "/parameter.txt");
  if (myfile.is_open()) {
    myfile << "Mem3DG Version: " << MEM3DG_VERSION << "\n";
    myfile << "Input Mesh:     " << inputMesh << "\n";
    myfile << "Physical parameters used: \n";
    myfile << "\n";
    myfile << "Kb:     " << system.parameters.bending.Kb << "\n"
           << "Kbc:   " << system.parameters.bending.Kbc << "\n"
           << "H0c:     " << system.parameters.bending.H0c << "\n"
           << "Kse:    " << system.meshProcessor.meshRegularizer.Kse << "\n"
           << "Ksl:    " << system.meshProcessor.meshRegularizer.Ksl << "\n"
           << "Kst:    " << system.meshProcessor.meshRegularizer.Kst << "\n"
           << "Ksg:    " << system.parameters.tension.Ksg << "\n"
           << "Kv:     " << system.parameters.osmotic.Kv << "\n"
           << "gamma:  " << system.parameters.dpd.gamma << "\n"
           << "Vt:     " << system.parameters.osmotic.Vt << "\n"
           << "kt:     " << system.parameters.temperature << "\n"
           << "Kf:     " << system.parameters.external.Kf << "\n";

    myfile << "\n";
    myfile << "Integration parameters used: \n";
    myfile << "\n";
    myfile << "dt:       " << timeStep << "\n"
           << "T:        " << totalTime << "\n"
           << "eps:		   " << tolerance << "\n"
           << "tSave:    " << savePeriod << "\n";
    myfile.close();

  } else
    std::cout << "Unable to open file";
}

void Integrator::getStatusLog(std::string nameOfFile, std::size_t frame,
                              double areaError, double volumeError,
                              double bendingError, double faceError,
                              std::string inputMesh) {
  std::ofstream myfile(nameOfFile);
  if (myfile.is_open()) {
    myfile << "Input Mesh: " << inputMesh << "\n";
    myfile << "Final parameter: \n";
    myfile << "\n";
    myfile << "Kb:     " << system.parameters.bending.Kb << "\n"
           << "Kbc:   " << system.parameters.bending.Kbc << "\n"
           << "H0c:     " << system.parameters.bending.H0c << "\n"
           << "Kse:    " << system.meshProcessor.meshRegularizer.Kse << "\n"
           << "Ksl:    " << system.meshProcessor.meshRegularizer.Ksl << "\n"
           << "Kst:    " << system.meshProcessor.meshRegularizer.Kst << "\n"
           << "Ksg:    " << system.parameters.tension.Ksg << "\n"
           << "Kv:     " << system.parameters.osmotic.Kv << "\n"
           << "gamma:  " << system.parameters.dpd.gamma << "\n"
           << "Vt:     " << system.parameters.osmotic.Vt << "\n"
           << "kt:     " << system.parameters.temperature << "\n"
           << "Kf:   " << system.parameters.external.Kf << "\n";

    myfile << "\n";
    myfile << "Integration: \n";
    myfile << "\n";
    myfile << "dt:    " << timeStep << "\n"
           << "T:     " << system.time << "\n"
           << "Frame: " << frame << "\n";

    myfile << "\n";
    myfile << "States: \n";
    myfile << "\n";
    myfile << "Bending Energy:   " << system.energy.bendingEnergy << "\n"
           << "Surface Energy:   " << system.energy.surfaceEnergy << "\n"
           << "Pressure Work:    " << system.energy.pressureEnergy << "\n"
           << "Kinetic Work:    " << system.energy.kineticEnergy << "\n"
           << "Adsorption Energy:  " << system.energy.adsorptionEnergy << "\n"
           << "Line tension Energy:  " << system.energy.dirichletEnergy << "\n"
           << "Total Energy:     " << system.energy.totalEnergy << "\n"
           << "Mech error norm:    " << system.mechErrorNorm << "\n"
           << "Chem error norm:    " << system.chemErrorNorm << "\n"
           << "\n"
           << "Surface area:     " << system.surfaceArea << " = "
           << system.surfaceArea / system.parameters.tension.At
           << " target surface area"
           << "\n"
           << "COM (x, y, z):		 "
           << toMatrix(system.vpg->inputVertexPositions).colwise().sum() /
                  system.vpg->inputVertexPositions.raw().rows()
           << "\n";

    myfile << "\n";
    myfile << "Errors: \n";
    myfile << "\n";
    myfile << "Bending error:       " << bendingError * 100 << "%"
           << "\n"
           << "Volume error:        " << volumeError * 100 << "%"
           << "\n"
           << "Surface area error:  " << areaError * 100 << "%"
           << "\n"
           << "Face area error:     " << faceError * 100 << "%"
           << "\n";

    myfile << "\n";
    myfile << "Options: \n";
    myfile << "\n";
    myfile << "Is considering protein: "
           << system.parameters.variation.isProteinVariation << "\n"
           << "Is vertex shift: "
           << system.meshProcessor.meshMutator.shiftVertex << "\n";

    myfile.close();
  } else
    std::cout << "Unable to open file";
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
