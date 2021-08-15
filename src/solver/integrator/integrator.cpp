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

double Integrator::backtrack(
    const double energy_pre,
    Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection, double rho,
    double c1) {

  auto physicalForceVec = toMatrix(f.forces.mechanicalForceVec);

  // validate the directions
  double positionProjection = 0;
  double chemicalProjection = 0;
  if (f.parameters.variation.isShapeVariation) {
    positionProjection =
        (physicalForceVec.array() * positionDirection.array()).sum();
    if (positionProjection < 0) {
      std::cout << "\nBacktracking line search: positional velocity on uphill "
                   "direction, use bare "
                   "gradient! \n"
                << std::endl;
      positionDirection = physicalForceVec;
      positionProjection =
          (physicalForceVec.array() * positionDirection.array()).sum();
    }
  }
  if (f.parameters.variation.isProteinVariation) {
    chemicalProjection =
        (f.forces.chemicalPotential.raw().array() * chemicalDirection.array())
            .sum();
    if (chemicalProjection < 0) {
      std::cout << "\nBacktracking line search: chemical direction on "
                   "uphill direction, "
                   "use bare "
                   "gradient! \n"
                << std::endl;
      chemicalDirection = f.forces.chemicalPotential.raw();
      chemicalProjection =
          (f.forces.chemicalPotential.raw().array() * chemicalDirection.array())
              .sum();
    }
  }

  // calculate initial energy as reference level
  const Eigen::Matrix<double, Eigen::Dynamic, 3> initial_pos =
      gc::EigenMap<double, 3>(f.vpg->inputVertexPositions);
  const Eigen::Matrix<double, Eigen::Dynamic, 1> initial_protein =
      f.proteinDensity.raw();
  const double init_time = f.time;

  // declare variables used in backtracking iterations
  double alpha = dt;
  std::size_t count = 0;

  // zeroth iteration
  if (f.parameters.variation.isShapeVariation) {
    toMatrix(f.vpg->inputVertexPositions) += alpha * positionDirection;
  }
  if (f.parameters.variation.isProteinVariation) {
    f.proteinDensity.raw() += alpha * chemicalDirection;
  }
  f.time += alpha;
  f.updateVertexPositions(false);
  f.computeFreeEnergy();

  while (true) {
    // Wolfe condition fulfillment
    if (f.energy.potE <
        (energy_pre - c1 * alpha * (positionProjection + chemicalProjection))) {
      break;
    }

    // limit of backtraking iterations
    if (alpha < 1e-5 * dt) {
      std::cout << "\nbacktrack: line search failure! Simulation "
                   "stopped. \n"
                << std::endl;
      lineSearchErrorBacktrack(alpha, initial_pos, initial_protein, true);
      EXIT = true;
      SUCCESS = false;
      break;
    }

    // backtracking time step
    alpha *= rho;
    if (f.parameters.variation.isShapeVariation) {
      toMatrix(f.vpg->inputVertexPositions) =
          initial_pos + alpha * positionDirection;
    }
    if (f.parameters.variation.isProteinVariation) {
      f.proteinDensity.raw() = initial_protein + alpha * chemicalDirection;
    }
    f.time = init_time + alpha;
    f.updateVertexPositions(false);
    f.computeFreeEnergy();

    // count the number of iterations
    count++;
  }

  // report the backtracking if verbose
  if (alpha != dt && verbosity > 3) {
    std::cout << "alpha: " << dt << " -> " << alpha << std::endl;
    std::cout << "mech norm: " << f.mechErrorNorm << std::endl;
    std::cout << "chem norm: " << f.chemErrorNorm << std::endl;
  }

  // If needed to test force-energy test
  const bool isDebug = false;
  if (isDebug) {
    lineSearchErrorBacktrack(alpha, initial_pos, initial_protein, isDebug);
  }

  return alpha;
}

void Integrator::lineSearchErrorBacktrack(
    const double &alpha, const EigenVectorX3dr &current_pos,
    const EigenVectorX1d &current_proteinDensity, bool runAll) {
  std::cout << "\nlineSearchErrorBacktracking ..." << std::endl;

  if (runAll || f.energy.potE > previousE.potE) {
    if (runAll || f.energy.BE > previousE.BE) {
      std::cout << "\nWith F_tol, BE has increased "
                << f.energy.BE - previousE.BE << " from " << previousE.BE
                << " to " << f.energy.BE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos +
          alpha * f.forces.maskForce(toMatrix(f.forces.bendingForceVec));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.BE > previousE.BE) {
        std::cout << "With only bending force, BE has increased "
                  << f.energy.BE - previousE.BE << " from " << previousE.BE
                  << " to " << f.energy.BE << ", expected dBE: "
                  << -alpha *
                         f.forces.maskForce(toMatrix(f.forces.bendingForceVec))
                             .squaredNorm()
                  << std::endl;
        // << -alpha *
        //        (f.F.mask(toMatrix(f.F.bendingForceVec)).array() *
        //         f.F.mask(toMatrix(f.F.bendingForceVec)).array())
        //            .sum()
        // << std::endl;
      }

      toMatrix(f.vpg->inputVertexPositions) = current_pos;
      f.proteinDensity.raw() =
          current_proteinDensity +
          alpha * f.parameters.Bc *
              f.forces.maskProtein(f.forces.bendingPotential.raw());
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.BE > previousE.BE) {
        std::cout << "With only bending potential, BE has increased "
                  << f.energy.BE - previousE.BE << " from " << previousE.BE
                  << " to " << f.energy.BE << ", expected dBE: "
                  << -alpha * f.parameters.Bc *
                         f.forces.maskProtein(f.forces.bendingPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }
    if (runAll || f.energy.sE > previousE.sE) {
      std::cout << "\nWith F_tol, sE has increased "
                << f.energy.sE - previousE.sE << " from " << previousE.sE
                << " to " << f.energy.sE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos +
          alpha * f.forces.maskForce(toMatrix(f.forces.capillaryForceVec));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.sE > previousE.sE) {
        std::cout << "With only capillary force, sE has increased "
                  << f.energy.sE - previousE.sE << " from " << previousE.sE
                  << " to " << f.energy.sE << ", expected dsE: "
                  << -alpha *
                         f.forces
                             .maskForce(toMatrix(f.forces.capillaryForceVec))
                             .squaredNorm()
                  << std::endl;
      }
    }
    if (runAll || f.energy.pE > previousE.pE) {
      std::cout << "\nWith F_tol, pE has increased "
                << f.energy.pE - previousE.pE << " from " << previousE.pE
                << " to " << f.energy.pE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos +
          alpha * f.forces.maskForce(toMatrix(f.forces.osmoticForceVec));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.pE > previousE.pE) {
        std::cout << "With only osmotic force, pE has increased "
                  << f.energy.pE - previousE.pE << " from " << previousE.pE
                  << " to " << f.energy.pE << ", expected dpE: "
                  << -alpha *
                         f.forces.maskForce(toMatrix(f.forces.osmoticForceVec))
                             .squaredNorm()
                  << std::endl;
      }
    }
    if (runAll || f.energy.aE > previousE.aE) {
      std::cout << "\nWith F_tol, aE has increased "
                << f.energy.aE - previousE.aE << " from " << previousE.aE
                << " to " << f.energy.aE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos +
          alpha * f.forces.maskForce(toMatrix(f.forces.adsorptionForceVec));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.aE > previousE.aE) {
        std::cout << "With only adsorption force, aE has increased "
                  << f.energy.aE - previousE.aE << " from " << previousE.aE
                  << " to " << f.energy.aE << ", expected daE: "
                  << -alpha *
                         f.forces
                             .maskForce(toMatrix(f.forces.adsorptionForceVec))
                             .squaredNorm()
                  << std::endl;
      }

      toMatrix(f.vpg->inputVertexPositions) = current_pos;
      f.proteinDensity.raw() =
          current_proteinDensity +
          alpha * f.parameters.Bc *
              f.forces.maskProtein(f.forces.adsorptionPotential.raw());
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.aE > previousE.aE) {
        std::cout << "With only adsorption potential, aE has increased "
                  << f.energy.aE - previousE.aE << " from " << previousE.aE
                  << " to " << f.energy.aE << ", expected dBE: "
                  << -alpha * f.parameters.Bc *
                         f.forces
                             .maskProtein(f.forces.adsorptionPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }
    if (runAll || f.energy.dE > previousE.dE) {
      std::cout << "\nWith F_tol, dE has increased "
                << f.energy.dE - previousE.dE << " from " << previousE.dE
                << " to " << f.energy.dE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos +
          alpha * f.forces.maskForce(toMatrix(f.forces.lineCapillaryForceVec));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.dE > previousE.dE) {
        std::cout << "With only line tension force, dE has increased "
                  << f.energy.dE - previousE.dE << " from " << previousE.dE
                  << " to " << f.energy.dE << ", expected ddE: "
                  << -alpha * f.forces
                                  .maskForce(
                                      toMatrix(f.forces.lineCapillaryForceVec))
                                  .squaredNorm()
                  << std::endl;
      }

      toMatrix(f.vpg->inputVertexPositions) = current_pos;
      f.proteinDensity.raw() =
          current_proteinDensity +
          alpha * f.parameters.Bc *
              f.forces.maskProtein(f.forces.diffusionPotential.raw());
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.dE > previousE.dE) {
        std::cout << "With only diffusion potential, dE has increased "
                  << f.energy.dE - previousE.dE << " from " << previousE.dE
                  << " to " << f.energy.dE << ", expected ddE: "
                  << -alpha * f.parameters.Bc *
                         f.forces.maskProtein(f.forces.diffusionPotential.raw())
                             .squaredNorm()
                  << std::endl;
      }
    }
    if (runAll || f.energy.exE > previousE.exE) {
      std::cout << "\nWith F_tol, exE has increased "
                << f.energy.exE - previousE.exE << " from " << previousE.exE
                << " to " << f.energy.exE << std::endl;

      f.proteinDensity.raw() = current_proteinDensity;
      toMatrix(f.vpg->inputVertexPositions) =
          current_pos + alpha * f.forces.maskForce(f.forces.addNormal(
                                    f.forces.externalForce.raw()));
      f.updateVertexPositions(false);
      f.computeFreeEnergy();
      if (runAll || f.energy.exE > previousE.exE) {
        std::cout << "With only external force, exE has increased "
                  << f.energy.exE - previousE.exE << " from " << previousE.exE
                  << " to " << f.energy.exE << ", expected dexE: "
                  << -alpha * f.forces
                                  .maskForce(f.forces.addNormal(
                                      f.forces.externalForce.raw()))
                                  .squaredNorm()
                  << std::endl;
      }
    }
  }
  if (runAll || f.energy.kE > previousE.kE) {
    std::cout << "\nWith F_tol, kE has increased " << f.energy.kE - previousE.kE
              << " from " << previousE.kE << " to " << f.energy.kE << std::endl;
  }
}

void Integrator::finitenessErrorBacktrack() {

  if (!std::isfinite(dt)) {
    EXIT = true;
    SUCCESS = false;
    std::cout << "time step is not finite!" << std::endl;
  }

  if (!std::isfinite(f.mechErrorNorm)) {
    EXIT = true;
    SUCCESS = false;

    if (!std::isfinite(toMatrix(f.velocity).norm())) {
      std::cout << "Velocity is not finite!" << std::endl;
    }

    if (!std::isfinite(toMatrix(f.forces.mechanicalForceVec).norm())) {
      if (!std::isfinite(toMatrix(f.forces.capillaryForceVec).norm())) {
        std::cout << "Capillary force is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.bendingForceVec).norm())) {
        std::cout << "Bending force is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.osmoticForceVec).norm())) {
        std::cout << "Osmotic force is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.lineCapillaryForceVec).norm())) {
        std::cout << "Line capillary force is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.externalForce).norm())) {
        std::cout << "External force is not finite!" << std::endl;
      }
    }
  }

  if (!std::isfinite(f.chemErrorNorm)) {
    EXIT = true;
    SUCCESS = false;

    if (!std::isfinite(toMatrix(f.proteinVelocity).norm())) {
      std::cout << "Protein velocity is not finite!" << std::endl;
    }

    if (!std::isfinite(toMatrix(f.forces.chemicalPotential).norm())) {
      if (!std::isfinite(toMatrix(f.forces.bendingPotential).norm())) {
        std::cout << "Bending Potential is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.interiorPenaltyPotential).norm())) {
        std::cout << "Protein interior penalty potential is not finite!"
                  << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.diffusionPotential).norm())) {
        std::cout << "Diffusion potential is not finite!" << std::endl;
      }
      if (!std::isfinite(toMatrix(f.forces.adsorptionPotential).norm())) {
        std::cout << "Adsorption potential is not finite!" << std::endl;
      }
    }
  }

  if (!std::isfinite(f.energy.totalE)) {
    EXIT = true;
    SUCCESS = false;
    if (!std::isfinite(f.energy.kE)) {
      std::cout << "Kinetic energy is not finite!" << std::endl;
    }
    if (!std::isfinite(f.energy.potE)) {
      if (!std::isfinite(f.energy.BE)) {
        std::cout << "Bending energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.sE)) {
        std::cout << "Surface energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.pE)) {
        std::cout << "Pressure energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.aE)) {
        std::cout << "Adsorption energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.dE)) {
        std::cout << "Line tension energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.exE)) {
        std::cout << "External force energy is not finite!" << std::endl;
      }
      if (!std::isfinite(f.energy.inE)) {
        std::cout << "Protein interior penalty energy is not finite!"
                  << std::endl;
      }
    }
  }
}

void Integrator::getForces() {
  auto vertexAngleNormal_e = gc::EigenMap<double, 3>(f.vpg->vertexNormals);

  f.computePhysicalForces();

  if ((f.parameters.dpd.gamma != 0) || (f.parameters.temp != 0)) {
    f.computeDPDForces(dt);
    DPDForce = rowwiseDotProduct(
        f.forces.maskForce(EigenMap<double, 3>(f.forces.dampingForce) +
                           EigenMap<double, 3>(f.forces.stochasticForce)),
        vertexAngleNormal_e);
  }

  // if (!f.mesh->hasBoundary()) {
  //   removeTranslation(physicalForceVec);
  //   removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
  //                  physicalForceVec);
  //   // removeTranslation(DPDPressure);
  //   // removeRotation(EigenMap<double, 3>(f.vpg->inputVertexPositions),
  //   // DPDPressure);
  // }
}

void Integrator::pressureConstraintThreshold(bool &EXIT,
                                             const bool isAugmentedLagrangian,
                                             const double dArea,
                                             const double ctol,
                                             double increment) {
  if (f.mechErrorNorm < tol && f.chemErrorNorm < tol) {
    if (isAugmentedLagrangian) { // augmented Lagrangian method
      if (dArea < ctol) {        // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG] = [" << f.parameters.lambdaSG << ", "
                  << "]";
        f.parameters.lambdaSG += f.parameters.tension.Ksg *
                                 (f.surfaceArea - f.refSurfaceArea) /
                                 f.refSurfaceArea;
        std::cout << " -> [" << f.parameters.lambdaSG << "]" << std::endl;
      }
    } else {              // incremental harmonic penalty method
      if (dArea < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[Ksg] = [" << f.parameters.tension.Ksg << "]";
        f.parameters.tension.Ksg *= increment;
        std::cout << " -> [" << f.parameters.tension.Ksg << "]" << std::endl;
      }
    }
  }
}

void Integrator::reducedVolumeThreshold(bool &EXIT,
                                        const bool isAugmentedLagrangian,
                                        const double dArea,
                                        const double dVolume, const double ctol,
                                        double increment) {
  if (f.mechErrorNorm < tol && f.chemErrorNorm < tol) {
    if (isAugmentedLagrangian) {            // augmented Lagrangian method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      } else { // iterate if not
        std::cout << "\n[lambdaSG, lambdaV] = [" << f.parameters.lambdaSG
                  << ", " << f.parameters.lambdaV << "]";
        f.parameters.lambdaSG += f.parameters.tension.Ksg *
                                 (f.surfaceArea - f.refSurfaceArea) /
                                 f.refSurfaceArea;
        f.parameters.lambdaV += f.parameters.osmotic.Kv *
                                (f.volume - f.parameters.osmotic.Vt) /
                                f.parameters.osmotic.Vt;
        std::cout << " -> [" << f.parameters.lambdaSG << ", "
                  << f.parameters.lambdaV << "]" << std::endl;
      }
    } else { // incremental harmonic penalty method
      if (dArea < ctol && dVolume < ctol) { // exit if fulfilled all constraints
        std::cout << "\nError norm smaller than tolerance." << std::endl;
        EXIT = true;
      }

      // iterate if not
      if (dArea > ctol) {
        std::cout << "\n[Ksg] = [" << f.parameters.tension.Ksg << "]";
        f.parameters.tension.Ksg *= 1.3;
        std::cout << " -> [" << f.parameters.tension.Ksg << "]" << std::endl;
      }
      if (dVolume > ctol) {
        std::cout << "\n[Kv] = [" << f.parameters.osmotic.Kv << "]";
        f.parameters.osmotic.Kv *= 1.3;
        std::cout << " -> [" << f.parameters.osmotic.Kv << "]" << std::endl;
      }
    }
  }
}

void Integrator::createNetcdfFile() {
  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  if (verbosity > 0) {
    fd.createNewFile(outputDir + "/" + trajFileName, *f.mesh, *f.refVpg,
                     TrajFile::NcFile::replace);
    fd.writeMask(toMatrix(f.forces.forceMask).rowwise().sum());
    if (!f.mesh->hasBoundary()) {
      fd.writeRefSurfArea(f.refSurfaceArea);
    }
  }
#endif
}

void Integrator::saveData() {
  // save variable to richData and save ply file
  if ((verbosity > 3 && !f.meshProcessor.meshMutator.isSplitEdge &&
       !f.meshProcessor.meshMutator.isCollapseEdge) ||
      (verbosity > 0 && (f.meshProcessor.meshMutator.isSplitEdge ||
                         f.meshProcessor.meshMutator.isCollapseEdge))) {
    char buffer[50];
    if (isJustGeometryPly) {
      sprintf(buffer, "/frame%d.obj", (int)frame);
    } else {
      sprintf(buffer, "/frame%d.ply", (int)frame);
    }
    f.saveRichData(outputDir + "/" + std::string(buffer), isJustGeometryPly);
  }

#ifdef MEM3DG_WITH_NETCDF
  // save variable to netcdf traj file
  if (verbosity > 0) {
    saveNetcdfData();
  }
#endif

  // print in-progress information in the console
  if (verbosity > 1) {
    std::cout << "\n"
              << "t: " << f.time << ", "
              << "n: " << frame << ", "
              << "isSmooth: " << f.isSmooth << "\n"
              << "dA/Area: " << dArea << "/" << f.surfaceArea << ", "
              << "dVP/Volume: " << dVP << "/" << f.volume << ", "
              << "h: "
              << toMatrix(f.vpg->inputVertexPositions).col(2).maxCoeff() << "\n"
              << "E_total: " << f.energy.totalE << "\n"
              << "E_pot: " << f.energy.potE << "\n"
              << "|e|Mech: " << f.mechErrorNorm << "\n"
              << "|e|Chem: " << f.chemErrorNorm << "\n"
              << "H: ["
              << (f.vpg->vertexMeanCurvatures.raw().array() /
                  f.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (f.vpg->vertexMeanCurvatures.raw().array() /
                  f.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]"
              << "\n"
              << "K: ["
              << (f.vpg->vertexGaussianCurvatures.raw().array() /
                  f.vpg->vertexDualAreas.raw().array())
                     .minCoeff()
              << ","
              << (f.vpg->vertexGaussianCurvatures.raw().array() /
                  f.vpg->vertexDualAreas.raw().array())
                     .maxCoeff()
              << "]" << std::endl;
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
                << ", and data saved to " + outputDir << std::endl;
      if (verbosity > 2) {
        f.saveRichData(outputDir + "/out.ply");
      }
    }
  }

  frame++;
}

void Integrator::markFileName(std::string marker_str) {
  std::string dirPath = outputDir;

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
void Integrator::saveNetcdfData() {
  std::size_t idx = fd.nFrames();

  // scalar quantities
  // write time
  fd.writeTime(idx, f.time);
  fd.writeIsSmooth(idx, f.isSmooth);
  // write geometry
  fd.writeVolume(idx, f.volume);
  fd.writeSurfArea(idx, f.mesh->hasBoundary() ? f.surfaceArea - f.refSurfaceArea
                                              : f.surfaceArea);
  fd.writeHeight(idx, toMatrix(f.vpg->inputVertexPositions).col(2).maxCoeff());
  // write energies
  fd.writeBendEnergy(idx, f.energy.BE);
  fd.writeSurfEnergy(idx, f.energy.sE);
  fd.writePressEnergy(idx, f.energy.pE);
  fd.writeKineEnergy(idx, f.energy.kE);
  fd.writeAdspEnergy(idx, f.energy.aE);
  fd.writeLineEnergy(idx, f.energy.dE);
  fd.writeTotalEnergy(idx, f.energy.totalE);
  // write Norms
  fd.writeErrorNorm(idx, f.mechErrorNorm);
  fd.writeChemErrorNorm(idx, f.chemErrorNorm);
  fd.writeBendNorm(idx, f.computeNorm(f.forces.bendingForce.raw()));
  fd.writeSurfNorm(idx, f.computeNorm(f.forces.capillaryForce.raw()));
  fd.writePressNorm(idx, f.computeNorm(f.forces.osmoticForce.raw()));
  fd.writeLineNorm(idx, f.computeNorm(f.forces.lineCapillaryForce.raw()));

  // vector quantities
  if (!f.meshProcessor.meshMutator.isSplitEdge &&
      !f.meshProcessor.meshMutator.isCollapseEdge) {
    // write velocity
    fd.writeVelocity(idx, EigenMap<double, 3>(f.velocity));
    // write protein density distribution
    fd.writeProteinDensity(idx, f.proteinDensity.raw());

    // write geometry
    fd.writeCoords(idx, EigenMap<double, 3>(f.vpg->inputVertexPositions));
    fd.writeTopoFrame(idx, f.mesh->getFaceVertexMatrix<std::uint32_t>());
    fd.writeMeanCurvature(idx, f.vpg->vertexMeanCurvatures.raw().array() /
                                   f.vpg->vertexDualAreas.raw().array());
    fd.writeGaussCurvature(idx, f.vpg->vertexGaussianCurvatures.raw().array() /
                                    f.vpg->vertexDualAreas.raw().array());
    fd.writeSponCurvature(idx, f.H0.raw());
    // fd.writeAngles(idx, f.vpg.cornerAngles.raw());
    // fd.writeH_H0_diff(idx,
    //                   ((f.H - f.H0).array() * (f.H -
    //                   f.H0).array()).matrix());

    // write pressures
    fd.writeBendingForce(idx, f.forces.bendingForce.raw());
    fd.writeCapillaryForce(idx, f.forces.capillaryForce.raw());
    fd.writeLineForce(idx, f.forces.lineCapillaryForce.raw());
    fd.writeOsmoticForce(idx, f.forces.osmoticForce.raw());
    fd.writeExternalForce(idx, f.forces.externalForce.raw());
    fd.writePhysicalForce(idx, f.forces.mechanicalForce.raw());
    fd.writeChemicalPotential(idx, f.forces.chemicalPotential.raw());
  }
}
#endif

void Integrator::getParameterLog(std::string inputMesh) {
  std::ofstream myfile(outputDir + "/parameter.txt");
  if (myfile.is_open()) {
    myfile << "Mem3DG Version: " << MEM3DG_VERSION << "\n";
    myfile << "Input Mesh:     " << inputMesh << "\n";
    myfile << "Physical parameters used: \n";
    myfile << "\n";
    myfile << "Kb:     " << f.parameters.bending.Kb << "\n"
           << "Kbc:   " << f.parameters.bending.Kbc << "\n"
           << "H0c:     " << f.parameters.bending.H0c << "\n"
           << "Kse:    " << f.meshProcessor.meshRegularizer.Kse << "\n"
           << "Ksl:    " << f.meshProcessor.meshRegularizer.Ksl << "\n"
           << "Kst:    " << f.meshProcessor.meshRegularizer.Kst << "\n"
           << "Ksg:    " << f.parameters.tension.Ksg << "\n"
           << "Kv:     " << f.parameters.osmotic.Kv << "\n"
           << "gamma:  " << f.parameters.dpd.gamma << "\n"
           << "Vt:     " << f.parameters.osmotic.Vt << "\n"
           << "kt:     " << f.parameters.temp << "\n"
           << "Kf:     " << f.parameters.external.Kf << "\n"
           << "conc:   " << f.parameters.external.conc << "\n"
           << "height: " << f.parameters.external.height << "\n";

    myfile << "\n";
    myfile << "Integration parameters used: \n";
    myfile << "\n";
    myfile << "dt:       " << dt << "\n"
           << "T:        " << total_time << "\n"
           << "eps:		   " << tol << "\n"
           << "tSave:    " << tSave << "\n";
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
    myfile << "Kb:     " << f.parameters.bending.Kb << "\n"
           << "Kbc:   " << f.parameters.bending.Kbc << "\n"
           << "H0c:     " << f.parameters.bending.H0c << "\n"
           << "Kse:    " << f.meshProcessor.meshRegularizer.Kse << "\n"
           << "Ksl:    " << f.meshProcessor.meshRegularizer.Ksl << "\n"
           << "Kst:    " << f.meshProcessor.meshRegularizer.Kst << "\n"
           << "Ksg:    " << f.parameters.tension.Ksg << "\n"
           << "Kv:     " << f.parameters.osmotic.Kv << "\n"
           << "gamma:  " << f.parameters.dpd.gamma << "\n"
           << "Vt:     " << f.parameters.osmotic.Vt << "\n"
           << "kt:     " << f.parameters.temp << "\n"
           << "Kf:   " << f.parameters.external.Kf << "\n"
           << "conc:   " << f.parameters.external.conc << "\n";

    myfile << "\n";
    myfile << "Integration: \n";
    myfile << "\n";
    myfile << "dt:    " << dt << "\n"
           << "T:     " << f.time << "\n"
           << "Frame: " << frame << "\n";

    myfile << "\n";
    myfile << "States: \n";
    myfile << "\n";
    myfile << "Bending Energy:   " << f.energy.BE << "\n"
           << "Surface Energy:   " << f.energy.sE << "\n"
           << "Pressure Work:    " << f.energy.pE << "\n"
           << "Kinetic Work:    " << f.energy.kE << "\n"
           << "Adsorption Energy:  " << f.energy.aE << "\n"
           << "Line tension Energy:  " << f.energy.dE << "\n"
           << "Total Energy:     " << f.energy.totalE << "\n"
           << "Mech error norm:    " << f.mechErrorNorm << "\n"
           << "Chem error norm:    " << f.chemErrorNorm << "\n"
           << "\n"
           << "Surface area:     " << f.surfaceArea << " = "
           << f.surfaceArea / f.refSurfaceArea << " target surface area"
           << "\n"
           << "COM (x, y, z):		 "
           << gc::EigenMap<double, 3>(f.vpg->inputVertexPositions)
                      .colwise()
                      .sum() /
                  f.vpg->inputVertexPositions.raw().rows()
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
           << f.parameters.variation.isProteinVariation << "\n"
           << "Is vertex shift: " << f.meshProcessor.meshMutator.shiftVertex << "\n";

    myfile.close();
  } else
    std::cout << "Unable to open file";
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
