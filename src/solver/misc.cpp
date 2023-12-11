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

#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::prescribeGeodesicMasks() {
  // Initialize the constant mask based on distance from the point specified
  if (parameters.variation.geodesicMask >
          geometry.geodesicDistance.raw().maxCoeff() ||
      parameters.variation.geodesicMask <
          geometry.geodesicDistance.raw().minCoeff()) {
    mem3dg_runtime_error("either all vertices or none are "
                         "in integration disk, "
                         "set radius = -1 to disable!");
  }
  for (gcs::Vertex v : geometry.mesh->vertices()) {
    forces.forceMask[v] =
        (geometry.geodesicDistance[v] < parameters.variation.geodesicMask)
            ? gc::Vector3{1, 1, 1}
            : gc::Vector3{0, 0, 0};
    forces.proteinMask[v] =
        (geometry.geodesicDistance[v] < parameters.variation.geodesicMask) ? 1
                                                                           : 0;
  }
};

void System::backtraceEnergyGrowth(const double timeStep,
                                   const Energy previousEnergy) {
  if (parameters.external.form != NULL)
    computeExternalWork(time, timeStep);
  computeTotalEnergy();
  std::cout << "<<<<<" << std::endl;
  mem3dg_runtime_warning("backtracing the component that leads to the energy "
                         "growth (<--'s highlight growing components)");
  std::cout << "\nUsing time step " << timeStep << ":" << std::endl;
  std::cout << "\n";
  auto checkGrowth = [](const double &previousEnergy,
                        const double &currentEnergy, std::string key) {
    if (currentEnergy != previousEnergy) {
      std::string marker = "";
      if (currentEnergy > previousEnergy || !std::isfinite(currentEnergy))
        marker = " <--";
      std::cout << key << " increased "
                << currentEnergy - previousEnergy
                // << " from " << previousEnergy << " to " << currentEnergy
                << marker << std::endl;
    }
  };
  checkGrowth(previousEnergy.potentialEnergy, energy.potentialEnergy,
              "potentialEnergy");
  checkGrowth(previousEnergy.spontaneousCurvatureEnergy,
              energy.spontaneousCurvatureEnergy, "spontaneousCurvatureEnergy");
  checkGrowth(previousEnergy.areaDifferenceEnergy, energy.areaDifferenceEnergy,
              "areaDifferenceEnergy");
  checkGrowth(previousEnergy.surfaceEnergy, energy.surfaceEnergy,
              "surfaceEnergy");
  checkGrowth(previousEnergy.pressureEnergy, energy.pressureEnergy,
              "pressureEnergy");
  checkGrowth(previousEnergy.adsorptionEnergy, energy.adsorptionEnergy,
              "adsorptionEnergy");
  checkGrowth(previousEnergy.aggregationEnergy, energy.aggregationEnergy,
              "aggregationEnergy");
  checkGrowth(previousEnergy.entropyEnergy, energy.entropyEnergy,
              "entropyEnergy");
  checkGrowth(previousEnergy.dirichletEnergy, energy.dirichletEnergy,
              "dirichletEnergy");
  checkGrowth(previousEnergy.externalWork, energy.externalWork, "externalWork");
  checkGrowth(previousEnergy.proteinInteriorPenalty,
              energy.proteinInteriorPenalty, "proteinInteriorPenalty");
  checkGrowth(previousEnergy.selfAvoidancePenalty, energy.selfAvoidancePenalty,
              "selfAvoidancePenalty");
  checkGrowth(previousEnergy.edgeSpringEnergy, energy.edgeSpringEnergy,
              "edgeSpringEnergy");
  checkGrowth(previousEnergy.faceSpringEnergy, energy.faceSpringEnergy,
              "faceSpringEnergy");
  checkGrowth(previousEnergy.lcrSpringEnergy, energy.lcrSpringEnergy,
              "lcrSpringEnergy");
  std::cout << ">>>>>" << std::endl;
};

bool System::testConservativeForcing(const double timeStep) {
  std::cout << "<<<<<" << std::endl;
  mem3dg_runtime_warning("numerically testing componentwise forcing--energy "
                         "relation (<--'s highlight possible inconsistency)");
  std::cout << "\nUsing time step " << timeStep << ":" << std::endl;
  std::cout << "\n";
  updateConfigurations();
  computeTotalEnergy();
  computeConservativeForcing();
  const Energy previousEnergy{energy};
  const EigenVectorX3dr previousPosition =
      toMatrix(geometry.vpg->inputVertexPositions);
  const EigenVectorX1d previousProteinDensity = proteinDensity.raw();

  bool SUCCESS = true;

  auto summerize = [](double actualDelta1, double expectedDelta1,
                      double actualDelta2, double expectedDelta2,
                      double stepFold, std::string key, double currentEnergy) {
    double expectRate = 2;
    std::array<std::size_t, 2> alignment{35, 60};
    std::array<std::size_t, 2> space;
    key = key + std::string(": ");
    space[0] = alignment[0] - key.length();

    std::ostringstream half1;
    half1 << "decrease " << actualDelta1 << ", ";
    space[1] = alignment[1] - half1.str().length() - alignment[0];

    std::string marker = " <-----\n";

    std::cout << key;
    for (std::size_t i = 0; i < space[0]; i++) // for console readability
      std::cout << " ";
    double difference_h = abs(expectedDelta1 - actualDelta1);
    if (actualDelta1 != 0) {
      std::cout << half1.str();
      for (std::size_t i = 0; i < space[1]; i++) // for console readability
        std::cout << " ";
      std::cout << "expect " << expectedDelta1;
      if ((actualDelta1 < 0) || !std::isfinite(currentEnergy) ||
          (difference_h / abs(actualDelta1)) > 0.01) { // criteria 1: closeness
        std::cout << marker << std::endl;
        return false;
      } else {
        std::cout << "" << std::endl;
      }
    } else {
      std::cout << "inactivated" << marker << std::endl;
      return false;
    }

    // test using x * timeStep for convergence
    space[0] = alignment[0];
    double difference_xh = abs(expectedDelta2 - actualDelta2);
    if (difference_h < 1e-15 && difference_xh < 1e-15) { // probably exact
      for (std::size_t i = 0; i < space[0]; i++) // for console readability
        std::cout << " ";
      std::cout << "exact\n" << std::endl;
      return true;
    }
    for (std::size_t i = 0; i < space[0]; i++) // for console readability
      std::cout << " ";
    std::ostringstream half2;
    half2 << "converge with "
          << pow(difference_xh / difference_h, 1 / expectRate) << ", ";
    space[1] = alignment[1] - half2.str().length() - alignment[0];
    std::cout << half2.str();
    for (std::size_t i = 0; i < space[1]; i++) // for console readability
      std::cout << " ";
    std::cout << "expect " << expectRate;
    if (abs(difference_xh / difference_h - pow(stepFold, expectRate)) >
        0.03) { // criteria 2: convergence rate
      std::cout << marker << std::endl;
      return false;
    } else {
      std::cout << "" << std::endl;
    }

    std::cout << "" << std::endl;
    return true;
  };

  // lambda function to test mechanical forces
  auto testMechanical = [previousProteinDensity, previousPosition, timeStep,
                         this, summerize](
                            gcs::VertexData<gc::Vector3> &forceVec,
                            const double &previousEnergy,
                            const double &currentEnergy, std::string forceKey,
                            std::string energyKey) {
    double stepFold = 2;
    std::string key = forceKey; // + std::string("--") + energyKey;

    // lambda function to compute energy delta
    auto computeDelta = [&](double dt) {
      proteinDensity.raw() = previousProteinDensity;
      toMatrix(geometry.vpg->inputVertexPositions) =
          previousPosition + dt * forces.maskForce(toMatrix(forceVec));
      updateConfigurations();
      computeTotalEnergy();
      double actualDelta =
          -currentEnergy +
          previousEnergy; // note that currentEnergy is passed by referenced,
                          // which is modified by the previous line
      double expectedDelta =
          dt * forces.maskForce(toMatrix(forceVec)).squaredNorm();
      return std::tuple<double, double>(actualDelta, expectedDelta);
    };

    double actualDelta1, expectedDelta1, actualDelta2, expectedDelta2;
    std::tie(actualDelta1, expectedDelta1) = computeDelta(timeStep);
    std::tie(actualDelta2, expectedDelta2) = computeDelta(stepFold * timeStep);
    return summerize(actualDelta1, expectedDelta1, actualDelta2, expectedDelta2,
                     stepFold, key, currentEnergy);
  };

  auto testChemical = [previousProteinDensity, previousPosition, timeStep, this,
                       summerize](gcs::VertexData<double> &potential,
                                  const double &previousEnergy,
                                  const double &currentEnergy,
                                  std::string potentialKey,
                                  std::string energyKey) {
    double stepFold = 2;
    std::string key = potentialKey; // + std::string("--") + energyKey;

    // lambda function to compute energy delta
    auto computeDelta = [&](double dt) {
      toMatrix(geometry.vpg->inputVertexPositions) = previousPosition;
      proteinDensity.raw() =
          previousProteinDensity +
          dt * parameters.proteinMobility * forces.maskProtein(potential.raw());
      updateConfigurations();
      computeTotalEnergy();
      double actualDelta = -currentEnergy + previousEnergy;
      double expectedDelta = dt * parameters.proteinMobility *
                             forces.maskProtein(potential.raw()).squaredNorm();
      return std::tuple<double, double>(actualDelta, expectedDelta);
    };
    double actualDelta1, expectedDelta1, actualDelta2, expectedDelta2;
    std::tie(actualDelta1, expectedDelta1) = computeDelta(timeStep);
    std::tie(actualDelta2, expectedDelta2) = computeDelta(stepFold * timeStep);
    return summerize(actualDelta1, expectedDelta1, actualDelta2, expectedDelta2,
                     stepFold, key, currentEnergy);
  };

  // ==========================================================
  // ================  Spontaneous curvature ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.spontaneousCurvatureForceVec,
                           previousEnergy.spontaneousCurvatureEnergy,
                           energy.spontaneousCurvatureEnergy,
                           "spontaneousCurvatureForceVec",
                           "spontaneousCurvatureEnergy") &&
            SUCCESS; // leaving && success to the last to ensure running all
                     // the tests
  SUCCESS = testChemical(forces.spontaneousCurvaturePotential,
                         previousEnergy.spontaneousCurvatureEnergy,
                         energy.spontaneousCurvatureEnergy,
                         "spontaneousCurvaturePotential",
                         "spontaneousCurvatureEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================    Area Difference     ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.areaDifferenceForceVec,
                           previousEnergy.areaDifferenceEnergy,
                           energy.areaDifferenceEnergy,
                           "areaDifferenceForceVec", "areaDifferenceEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================ Deviatoric curvature   ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.deviatoricCurvatureForceVec,
                           previousEnergy.deviatoricCurvatureEnergy,
                           energy.deviatoricCurvatureEnergy,
                           "deviatoricCurvatureForceVec",
                           "deviatoricCurvatureEnergy") &&
            SUCCESS;
  SUCCESS = testChemical(forces.deviatoricCurvaturePotential,
                         previousEnergy.deviatoricCurvatureEnergy,
                         energy.deviatoricCurvatureEnergy,
                         "deviatoricCurvaturePotential",
                         "deviatoricCurvatureEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================      Capillary         ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.capillaryForceVec,
                           previousEnergy.surfaceEnergy, energy.surfaceEnergy,
                           "capillaryForceVec", "surfaceEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================      Osmotic           ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.osmoticForceVec,
                           previousEnergy.pressureEnergy, energy.pressureEnergy,
                           "osmoticForceVec", "pressureEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================      Adsorption        ==================
  // ==========================================================
  SUCCESS =
      testMechanical(forces.adsorptionForceVec, previousEnergy.adsorptionEnergy,
                     energy.adsorptionEnergy, "adsorptionForceVec",
                     "adsorptionEnergy") &&
      SUCCESS;
  SUCCESS =
      testChemical(forces.adsorptionPotential, previousEnergy.adsorptionEnergy,
                   energy.adsorptionEnergy, "adsorptionPotential",
                   "adsorptionEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================      Entropy           ==================
  // ==========================================================
  SUCCESS =

      testMechanical(forces.entropyForceVec, previousEnergy.entropyEnergy,
                     energy.entropyEnergy, "entropyForceVec",
                     "entropyEnergy") &&
      SUCCESS;
  SUCCESS =

      testChemical(forces.entropyPotential, previousEnergy.entropyEnergy,
                   energy.entropyEnergy, "entropyPotential", "entropyEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================      Aggregation       ==================
  // ==========================================================
  SUCCESS =

      testMechanical(forces.aggregationForceVec,
                     previousEnergy.aggregationEnergy, energy.aggregationEnergy,
                     "aggregationForceVec", "aggregationEnergy") &&
      SUCCESS;
  SUCCESS =

      testChemical(forces.aggregationPotential,
                   previousEnergy.aggregationEnergy, energy.aggregationEnergy,
                   "aggregationPotential", "aggregationEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================      Dirichlet         ==================
  // ==========================================================
  SUCCESS =

      testMechanical(forces.lineCapillaryForceVec,
                     previousEnergy.dirichletEnergy, energy.dirichletEnergy,
                     "lineCapillaryForceVec", "dirichletEnergy") &&
      SUCCESS;
  SUCCESS = testChemical(forces.dirichletPotential,
                         previousEnergy.dirichletEnergy, energy.dirichletEnergy,
                         "dirichletPotential", "dirichletEnergy") &&
            SUCCESS;
  // ==========================================================
  // ================   Self-avoidance       ==================
  // ==========================================================
  SUCCESS = testMechanical(forces.selfAvoidanceForceVec,
                           previousEnergy.selfAvoidancePenalty,
                           energy.selfAvoidancePenalty, "selfAvoidanceForceVec",
                           "selfAvoidancePenalty") &&
            SUCCESS;
  // ==========================================================
  // ================        Edge spring     ==================
  // ==========================================================
  SUCCESS =
      testMechanical(forces.edgeSpringForceVec, previousEnergy.edgeSpringEnergy,
                     energy.edgeSpringEnergy, "edgeSpringForceVec",
                     "edgeSpringEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================        Face spring     ==================
  // ==========================================================
  SUCCESS =
      testMechanical(forces.faceSpringForceVec, previousEnergy.faceSpringEnergy,
                     energy.faceSpringEnergy, "faceSpringForceVec",
                     "faceSpringEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================        LCR spring      ==================
  // ==========================================================
  SUCCESS =
      testMechanical(forces.lcrSpringForceVec, previousEnergy.lcrSpringEnergy,
                     energy.lcrSpringEnergy, "lcrSpringForceVec",
                     "lcrSpringEnergy") &&
      SUCCESS;
  // ==========================================================
  // ================     Interior Penalty   ==================
  // ==========================================================
  SUCCESS =
      testChemical(forces.interiorPenaltyPotential,
                   previousEnergy.proteinInteriorPenalty,
                   energy.proteinInteriorPenalty, "interiorPenaltyPotential",
                   "proteinInteriorPenalty") &&
      SUCCESS;

  // recover the configuration
  toMatrix(geometry.vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() = previousProteinDensity;
  updateConfigurations();
  computeTotalEnergy();
  computeConservativeForcing();

  std::cout << ">>>>>" << std::endl;
  return SUCCESS;
}

void System::check_pcg() {
  // Generate a normal distribution around that mean
  std::normal_distribution<> normal_dist(0, 2);

  // Make a copy of the RNG state to use later
  pcg32 rng_checkpoint = rng;

  // Produce histogram
  std::map<int, int> hist;
  for (int n = 0; n < 10000; ++n) {
    ++hist[std::round(normal_dist(rng))];
  }
  std::cout << "Normal distribution around " << 0 << ":\n";
  for (auto p : hist) {
    std::cout << std::fixed << std::setprecision(1) << std::setw(2) << p.first
              << ' ' << std::string(p.second / 30, '*') << '\n';
  }

  // Produce information about RNG usage
  std::cout << "Required " << (rng - rng_checkpoint) << " random numbers.\n";
}

bool System::checkFiniteness() {
  bool finite = true;
  if (!std::isfinite(mechErrorNorm)) {
    finite = false;
    if (!std::isfinite(toMatrix(velocity).norm())) {
      mem3dg_runtime_warning("Velocity is not finite!");
    }

    if (!std::isfinite(toMatrix(forces.mechanicalForceVec).norm())) {
      if (!std::isfinite(toMatrix(forces.capillaryForceVec).norm())) {
        mem3dg_runtime_warning("Capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.adsorptionForceVec).norm())) {
        mem3dg_runtime_warning("Adsorption force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.aggregationForceVec).norm())) {
        mem3dg_runtime_warning("Aggregation force is not finite!");
      }
      if (!std::isfinite(
              toMatrix(forces.spontaneousCurvatureForceVec).norm())) {
        mem3dg_runtime_warning("Spontaneous curvature force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.deviatoricCurvatureForceVec).norm())) {
        mem3dg_runtime_warning("Deviatoric curvature force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.areaDifferenceForceVec).norm())) {
        mem3dg_runtime_warning("Area difference force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.osmoticForceVec).norm())) {
        mem3dg_runtime_warning("Osmotic force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.lineCapillaryForceVec).norm())) {
        mem3dg_runtime_warning("Line capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.externalForceVec).norm())) {
        mem3dg_runtime_warning("External force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.selfAvoidanceForceVec).norm())) {
        mem3dg_runtime_warning("Self avoidance force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.entropyForceVec).norm())) {
        mem3dg_runtime_warning("Entropy force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.springForceVec).norm())) {
        mem3dg_runtime_warning("Spring force is not finite!");
      }
    }
  }

  if (!std::isfinite(chemErrorNorm)) {
    finite = false;

    if (!std::isfinite(proteinRateOfChange.raw().norm())) {
      mem3dg_runtime_warning("Protein velocity is not finite!");
    }

    if (!std::isfinite(forces.chemicalPotential.raw().norm())) {
      if (!std::isfinite(forces.spontaneousCurvaturePotential.raw().norm())) {
        mem3dg_runtime_warning(
            "Spontaneous curvature Potential is not finite!");
      }
      if (!std::isfinite(forces.deviatoricCurvaturePotential.raw().norm())) {
        mem3dg_runtime_warning("Deviatoric curvature Potential is not finite!");
      }
      if (!std::isfinite(forces.interiorPenaltyPotential.raw().norm())) {
        mem3dg_runtime_warning(
            "Protein interior penalty potential is not finite!");
      }
      if (!std::isfinite(forces.dirichletPotential.raw().norm())) {
        mem3dg_runtime_warning("Dirichlet potential is not finite!");
      }
      if (!std::isfinite(forces.adsorptionPotential.raw().norm())) {
        mem3dg_runtime_warning("Adsorption potential is not finite!");
      }
      if (!std::isfinite(forces.aggregationPotential.raw().norm())) {
        mem3dg_runtime_warning("Aggregation potential is not finite!");
      }
      if (!std::isfinite(forces.entropyPotential.raw().norm())) {
        mem3dg_runtime_warning("Entropy potential is not finite!");
      }
    }
  }

  if (!std::isfinite(energy.totalEnergy)) {
    finite = false;
    if (!std::isfinite(energy.kineticEnergy)) {
      mem3dg_runtime_warning("Kinetic energy is not finite!");
    }
    if (!std::isfinite(energy.externalWork)) {
      mem3dg_runtime_warning("External work is not finite!");
    }
    if (!std::isfinite(energy.potentialEnergy)) {
      if (!std::isfinite(energy.spontaneousCurvatureEnergy)) {
        mem3dg_runtime_warning("Spontaneous curvature energy is not finite!");
      }
      if (!std::isfinite(energy.deviatoricCurvatureEnergy)) {
        mem3dg_runtime_warning("Deviatoric curvature energy is not finite!");
      }
      if (!std::isfinite(energy.areaDifferenceEnergy)) {
        mem3dg_runtime_warning("Area difference energy is not finite!");
      }
      if (!std::isfinite(energy.surfaceEnergy)) {
        mem3dg_runtime_warning("Surface energy is not finite!");
      }
      if (!std::isfinite(energy.pressureEnergy)) {
        mem3dg_runtime_warning("Pressure energy is not finite!");
      }
      if (!std::isfinite(energy.adsorptionEnergy)) {
        mem3dg_runtime_warning("Adsorption energy is not finite!");
      }
      if (!std::isfinite(energy.aggregationEnergy)) {
        mem3dg_runtime_warning("Aggregation energy is not finite!");
      }
      if (!std::isfinite(energy.dirichletEnergy)) {
        mem3dg_runtime_warning("Line tension energy is not finite!");
      }
      if (!std::isfinite(energy.proteinInteriorPenalty)) {
        mem3dg_runtime_warning(
            "Protein interior penalty energy is not finite!");
      }
      if (!std::isfinite(energy.selfAvoidancePenalty)) {
        mem3dg_runtime_warning(
            "Membrane self-avoidance penalty energy is not finite!");
      }
      if (!std::isfinite(energy.entropyEnergy)) {
        mem3dg_runtime_warning("Entropy energy is not finite!");
      }
      if (!std::isfinite(energy.lcrSpringEnergy)) {
        mem3dg_runtime_warning("lcr spring energy is not finite!");
      }
      if (!std::isfinite(energy.edgeSpringEnergy)) {
        mem3dg_runtime_warning("edge spring energy is not finite!");
      }
      if (!std::isfinite(energy.faceSpringEnergy)) {
        mem3dg_runtime_warning("face spring energy is not finite!");
      }
    }
  }

  return finite;
}

} // namespace solver
} // namespace mem3dg
