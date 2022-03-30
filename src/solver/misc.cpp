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

void System::testForceComputation(const double timeStep) {
  computeTotalEnergy();
  computePhysicalForcing();
  const Energy previousE{energy};
  const EigenVectorX3dr previousR = toMatrix(vpg->inputVertexPositions);
  const EigenVectorX1d previousPhi = toMatrix(proteinDensity);

  vpg->inputVertexPositions += forces.mechanicalForceVec * timeStep;
  proteinDensity += parameters.proteinMobility * forces.chemicalPotential /
                    vpg->vertexDualAreas * timeStep;

  updateConfigurations();
  computeTotalEnergy();
  testForceComputation(timeStep, previousR, previousPhi, previousE);

  toMatrix(vpg->inputVertexPositions) = previousR;
  toMatrix(proteinDensity) = previousPhi;
  updateConfigurations();
  computeTotalEnergy();
  computePhysicalForcing();
}

void System::updateGeodesicsDistance() {
  // update geodesic distance
  gcs::HeatMethodDistanceSolver heatSolver(*vpg);
  geodesicDistance = heatSolver.computeDistance(center);
}

void System::prescribeGeodesicProteinDensityDistribution() {
  std::array<double, 2> r_heter{
      parameters.protein.geodesicProteinDensityDistribution[0],
      parameters.protein.geodesicProteinDensityDistribution[1]};
  vpg->requireVertexTangentBasis();
  if (parameters.protein.profile == "gaussian") {
    gaussianDistribution(proteinDensity.raw(), geodesicDistance.raw(),
                         vpg->inputVertexPositions -
                             vpg->inputVertexPositions[center.nearestVertex()],
                         vpg->vertexTangentBasis[center.nearestVertex()],
                         r_heter);
  } else if (parameters.protein.profile == "tanh") {
    tanhDistribution(proteinDensity.raw(), geodesicDistance.raw(),
                     vpg->inputVertexPositions -
                         vpg->inputVertexPositions[center.nearestVertex()],
                     vpg->vertexTangentBasis[center.nearestVertex()],
                     parameters.protein.tanhSharpness, r_heter);
  }
  vpg->unrequireVertexTangentBasis();
  proteinDensity.raw() *=
      parameters.protein.geodesicProteinDensityDistribution[2] -
      parameters.protein.geodesicProteinDensityDistribution[3];
  proteinDensity.raw().array() +=
      parameters.protein.geodesicProteinDensityDistribution[3];
}

void System::prescribeGeodesicMasks() {
  // Initialize the constant mask based on distance from the point specified
  if (parameters.variation.geodesicMask != -1) {
    if (parameters.variation.geodesicMask > geodesicDistance.raw().maxCoeff() ||
        parameters.variation.geodesicMask < geodesicDistance.raw().minCoeff()) {
      mem3dg_runtime_error("either all vertices or none is "
                           "initializeConstantsin integration disk, "
                           "set radius = -1 to disable!");
    }
    for (gcs::Vertex v : mesh->vertices()) {
      forces.forceMask[v] =
          (geodesicDistance[v] < parameters.variation.geodesicMask)
              ? gc::Vector3{1, 1, 1}
              : gc::Vector3{0, 0, 0};
      forces.proteinMask[v] =
          (geodesicDistance[v] < parameters.variation.geodesicMask) ? 1 : 0;
    }
  }
};

void System::testForceComputation(const double timeStep,
                                  const EigenVectorX3dr previousPosition,
                                  const EigenVectorX1d previousProteinDensity,
                                  const Energy previousEnergy) {

  // cache the energy when applied the total force
  // proteinDensity.raw() = previousProteinDensity;
  // toMatrix(vpg->inputVertexPositions) =
  //     previousPosition +
  //     timeStep * forces.maskForce(
  //                 toMatrix(forces.osmoticForceVec) +
  //                 toMatrix(forces.capillaryForceVec) +
  //                 toMatrix(forces.bendingForceVec) +
  //                 toMatrix(forces.lineCapillaryForceVec) +
  //                 toMatrix(forces.adsorptionForceVec) +
  //                 toMatrix(forces.aggregationForceVec) +
  //                 toMatrix(forces.externalForceVec) +
  //                 toMatrix(forces.selfAvoidanceForceVec));
  // // * toMatrix(forces.mechanicalForceVec);
  // updateConfigurations();
  if (parameters.external.Kf != 0)
    computeExternalWork(time, timeStep);
  computeTotalEnergy();
  const Energy totalForceEnergy{energy};

  // ==========================================================
  // ================    Potential Energy    ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.potentialEnergy > previousEnergy.potentialEnergy)
    std::cout << "*";
  std::cout << "With F_tol, potE has increased "
            << totalForceEnergy.potentialEnergy - previousEnergy.potentialEnergy
            << " from " << previousEnergy.potentialEnergy << " to "
            << totalForceEnergy.potentialEnergy << std::endl;

  // ==========================================================
  // ================        Bending         ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.bendingEnergy > previousEnergy.bendingEnergy)
    std::cout << "*";
  std::cout << "With F_tol, BE has increased "
            << totalForceEnergy.bendingEnergy - previousEnergy.bendingEnergy
            << " from " << previousEnergy.bendingEnergy << " to "
            << totalForceEnergy.bendingEnergy << std::endl;

  // =======================================
  // =======      Mechanical        ========
  // =======================================
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.bendingForceVec));
  updateConfigurations();

  // test if bending energy increases
  computeBendingEnergy();
  if (energy.bendingEnergy > previousEnergy.bendingEnergy)
    std::cout << "*";
  std::cout << "With only bending force, BE has increased "
            << energy.bendingEnergy - previousEnergy.bendingEnergy << " from "
            << previousEnergy.bendingEnergy << " to " << energy.bendingEnergy
            << ", expected dBE: "
            << -timeStep * forces.maskForce(toMatrix(forces.bendingForceVec))
                               .squaredNorm()
            << std::endl;

  // =======================================
  // =======       Chemical         ========
  // =======================================
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() = previousProteinDensity +
                         timeStep * parameters.proteinMobility *
                             forces.maskProtein(forces.bendingPotential.raw());
  updateConfigurations();

  // test if bending energy increases
  computeBendingEnergy();
  if (energy.bendingEnergy > previousEnergy.bendingEnergy)
    std::cout << "*";
  std::cout
      << "With only bending potential, BE has increased "
      << energy.bendingEnergy - previousEnergy.bendingEnergy << " from "
      << previousEnergy.bendingEnergy << " to " << energy.bendingEnergy
      << ", expected dBE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.bendingPotential.raw()).squaredNorm()
      << std::endl;

  // ==========================================================
  // ================      Deviatoric        ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.deviatoricEnergy > previousEnergy.deviatoricEnergy)
    std::cout << "*";
  std::cout << "With F_tol, DevE has increased "
            << totalForceEnergy.deviatoricEnergy -
                   previousEnergy.deviatoricEnergy
            << " from " << previousEnergy.deviatoricEnergy << " to "
            << totalForceEnergy.deviatoricEnergy << std::endl;

  // =======================================
  // =======      Mechanical        ========
  // =======================================
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.deviatoricForceVec));
  updateConfigurations();

  // test if deviatoric energy increases
  computeDeviatoricEnergy();
  if (energy.deviatoricEnergy > previousEnergy.deviatoricEnergy)
    std::cout << "*";
  std::cout << "With only deviatoric force, DevE has increased "
            << energy.deviatoricEnergy - previousEnergy.deviatoricEnergy
            << " from " << previousEnergy.deviatoricEnergy << " to "
            << energy.deviatoricEnergy << ", expected dDevE: "
            << -timeStep * forces.maskForce(toMatrix(forces.deviatoricForceVec))
                               .squaredNorm()
            << std::endl;

  // =======================================
  // =======       Chemical         ========
  // =======================================
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() =
      previousProteinDensity +
      timeStep * parameters.proteinMobility *
          forces.maskProtein(forces.deviatoricPotential.raw());
  updateConfigurations();

  // test if deviatoric energy increases
  computeDeviatoricEnergy();
  if (energy.deviatoricEnergy > previousEnergy.deviatoricEnergy)
    std::cout << "*";
  std::cout
      << "With only deviatoric potential, DevE has increased "
      << energy.deviatoricEnergy - previousEnergy.deviatoricEnergy << " from "
      << previousEnergy.deviatoricEnergy << " to " << energy.deviatoricEnergy
      << ", expected dDevE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.deviatoricPotential.raw()).squaredNorm()
      << std::endl;

  // ==========================================================
  // ================      Capillary         ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.surfaceEnergy > previousEnergy.surfaceEnergy)
    std::cout << "*";
  std::cout << "With F_tol, sE has increased "
            << totalForceEnergy.surfaceEnergy - previousEnergy.surfaceEnergy
            << " from " << previousEnergy.surfaceEnergy << " to "
            << totalForceEnergy.surfaceEnergy << std::endl;

  // =======================================
  // =======      Mechanical        ========
  // =======================================
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.capillaryForceVec));
  updateConfigurations();
  computeSurfaceEnergy();
  if (energy.surfaceEnergy > previousEnergy.surfaceEnergy)
    std::cout << "*";
  std::cout << "With only capillary force, sE has increased "
            << energy.surfaceEnergy - previousEnergy.surfaceEnergy << " from "
            << previousEnergy.surfaceEnergy << " to " << energy.surfaceEnergy
            << ", expected dsE: "
            << -timeStep * forces.maskForce(toMatrix(forces.capillaryForceVec))
                               .squaredNorm()
            << std::endl;

  // ==========================================================
  // ================      Osmotic           ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.pressureEnergy > previousEnergy.pressureEnergy)
    std::cout << "*";
  std::cout << "With F_tol, pE has increased "
            << totalForceEnergy.pressureEnergy - previousEnergy.pressureEnergy
            << " from " << previousEnergy.pressureEnergy << " to "
            << totalForceEnergy.pressureEnergy << std::endl;

  // =======================================
  // =======      Mechanical        ========
  // =======================================
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.osmoticForceVec));
  updateConfigurations();
  computePressureEnergy();
  if (energy.pressureEnergy > previousEnergy.pressureEnergy)
    std::cout << "*";
  std::cout << "With only osmotic force, pE has increased "
            << energy.pressureEnergy - previousEnergy.pressureEnergy << " from "
            << previousEnergy.pressureEnergy << " to " << energy.pressureEnergy
            << ", expected dpE: "
            << -timeStep * forces.maskForce(toMatrix(forces.osmoticForceVec))
                               .squaredNorm()
            << std::endl;

  // ==========================================================
  // ================      Adsorption        ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.adsorptionEnergy > previousEnergy.adsorptionEnergy)
    std::cout << "*";
  std::cout << "With F_tol, aE has increased "
            << totalForceEnergy.adsorptionEnergy -
                   previousEnergy.adsorptionEnergy
            << " from " << previousEnergy.adsorptionEnergy << " to "
            << totalForceEnergy.adsorptionEnergy << std::endl;

  // =======================================
  // =======      Mechanical        ========
  // =======================================
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.adsorptionForceVec));
  updateConfigurations();
  computeAdsorptionEnergy();
  if (energy.adsorptionEnergy > previousEnergy.adsorptionEnergy)
    std::cout << "*";
  std::cout << "With only adsorption force, aE has increased "
            << energy.adsorptionEnergy - previousEnergy.adsorptionEnergy
            << " from " << previousEnergy.adsorptionEnergy << " to "
            << energy.adsorptionEnergy << ", expected daE: "
            << -timeStep * forces.maskForce(toMatrix(forces.adsorptionForceVec))
                               .squaredNorm()
            << std::endl;

  // =======================================
  // =======       Chemical         ========
  // =======================================
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() =
      previousProteinDensity +
      timeStep * parameters.proteinMobility *
          forces.maskProtein(forces.adsorptionPotential.raw());
  updateConfigurations();
  computeAdsorptionEnergy();
  if (energy.adsorptionEnergy > previousEnergy.adsorptionEnergy)
    std::cout << "*";
  std::cout
      << "With only adsorption potential, aE has increased "
      << energy.adsorptionEnergy - previousEnergy.adsorptionEnergy << " from "
      << previousEnergy.adsorptionEnergy << " to " << energy.adsorptionEnergy
      << ", expected daE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.adsorptionPotential.raw()).squaredNorm()
      << std::endl;

  // ==========================================================
  // ================      Aggregation       ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.aggregationEnergy > previousEnergy.aggregationEnergy)
    std::cout << "*";
  std::cout << "With F_tol, aggE has increased "
            << totalForceEnergy.aggregationEnergy -
                   previousEnergy.aggregationEnergy
            << " from " << previousEnergy.aggregationEnergy << " to "
            << totalForceEnergy.aggregationEnergy << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.aggregationForceVec));
  updateConfigurations();
  computeAggregationEnergy();
  if (energy.aggregationEnergy > previousEnergy.aggregationEnergy)
    std::cout << "*";
  std::cout << "With only aggregation force, aggE has increased "
            << energy.aggregationEnergy - previousEnergy.aggregationEnergy
            << " from " << previousEnergy.aggregationEnergy << " to "
            << energy.aggregationEnergy << ", expected daggE: "
            << -timeStep *
                   forces.maskForce(toMatrix(forces.aggregationForceVec))
                       .squaredNorm()
            << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() =
      previousProteinDensity +
      timeStep * parameters.proteinMobility *
          forces.maskProtein(forces.aggregationPotential.raw());
  updateConfigurations();
  computeAggregationEnergy();
  if (energy.aggregationEnergy > previousEnergy.aggregationEnergy)
    std::cout << "*";
  std::cout
      << "With only aggregation potential, aggE has increased "
      << energy.aggregationEnergy - previousEnergy.aggregationEnergy << " from "
      << previousEnergy.aggregationEnergy << " to " << energy.aggregationEnergy
      << ", expected daggE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.aggregationPotential.raw()).squaredNorm()
      << std::endl;

  // ==========================================================
  // ================      Dirichlet         ==================
  // ==========================================================
  std::cout << "\n";
  if (energy.dirichletEnergy > previousEnergy.dirichletEnergy)
    std::cout << "*";
  std::cout << "With F_tol, dE has increased "
            << energy.dirichletEnergy - previousEnergy.dirichletEnergy
            << " from " << previousEnergy.dirichletEnergy << " to "
            << energy.dirichletEnergy << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.lineCapillaryForceVec));
  updateConfigurations();
  computeDirichletEnergy();
  if (energy.dirichletEnergy > previousEnergy.dirichletEnergy)
    std::cout << "*";
  std::cout << "With only line tension force, dE has increased "
            << energy.dirichletEnergy - previousEnergy.dirichletEnergy
            << " from " << previousEnergy.dirichletEnergy << " to "
            << energy.dirichletEnergy << ", expected ddE: "
            << -timeStep *
                   forces.maskForce(toMatrix(forces.lineCapillaryForceVec))
                       .squaredNorm()
            << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() =
      previousProteinDensity +
      timeStep * parameters.proteinMobility *
          forces.maskProtein(forces.diffusionPotential.raw());
  updateConfigurations();
  computeDirichletEnergy();
  if (energy.dirichletEnergy > previousEnergy.dirichletEnergy)
    std::cout << "*";
  std::cout
      << "With only diffusion potential, dE has increased "
      << energy.dirichletEnergy - previousEnergy.dirichletEnergy << " from "
      << previousEnergy.dirichletEnergy << " to " << energy.dirichletEnergy
      << ", expected ddE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.diffusionPotential.raw()).squaredNorm()
      << std::endl;

  // ==========================================================
  // ================   Self-avoidance       ==================
  // ==========================================================
  std::cout << "\n";
  if (energy.selfAvoidancePenalty > previousEnergy.selfAvoidancePenalty)
    std::cout << "*";
  std::cout << "With F_tol, selfE has increased "
            << energy.selfAvoidancePenalty - previousEnergy.selfAvoidancePenalty
            << " from " << previousEnergy.selfAvoidancePenalty << " to "
            << energy.selfAvoidancePenalty << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  proteinDensity.raw() = previousProteinDensity;
  toMatrix(vpg->inputVertexPositions) =
      previousPosition +
      timeStep * forces.maskForce(toMatrix(forces.selfAvoidanceForceVec));
  updateConfigurations();
  computeSelfAvoidanceEnergy();
  if (energy.selfAvoidancePenalty > previousEnergy.selfAvoidancePenalty)
    std::cout << "*";
  std::cout << "With only self avoidance penalty force, selfE has increased "
            << energy.selfAvoidancePenalty - previousEnergy.selfAvoidancePenalty
            << " from " << previousEnergy.selfAvoidancePenalty << " to "
            << energy.selfAvoidancePenalty << ", expected dselfE: "
            << -timeStep *
                   forces.maskForce(toMatrix(forces.selfAvoidanceForceVec))
                       .squaredNorm()
            << std::endl;

  // ==========================================================
  // ================     Interior Penalty   ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.proteinInteriorPenalty >
      previousEnergy.proteinInteriorPenalty)
    std::cout << "*";
  std::cout << "With F_tol, inPE has increased "
            << totalForceEnergy.proteinInteriorPenalty -
                   previousEnergy.proteinInteriorPenalty
            << " from " << previousEnergy.proteinInteriorPenalty << " to "
            << totalForceEnergy.proteinInteriorPenalty << std::endl;

  // test single-force-energy computation
  // perturb the configuration
  toMatrix(vpg->inputVertexPositions) = previousPosition;
  proteinDensity.raw() =
      previousProteinDensity +
      timeStep * parameters.proteinMobility *
          forces.maskProtein(forces.interiorPenaltyPotential.raw());
  updateConfigurations();
  computeProteinInteriorPenalty();
  if (energy.proteinInteriorPenalty > previousEnergy.proteinInteriorPenalty)
    std::cout << "*";
  std::cout
      << "With only protein interior penalty potential, inPE has increased "
      << energy.proteinInteriorPenalty - previousEnergy.proteinInteriorPenalty
      << " from " << previousEnergy.proteinInteriorPenalty << " to "
      << energy.proteinInteriorPenalty << ", expected dinPE: "
      << -timeStep * parameters.proteinMobility *
             forces.maskProtein(forces.interiorPenaltyPotential.raw())
                 .squaredNorm()
      << std::endl;

  // ==========================================================
  // ================     External Force     ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.externalWork < previousEnergy.externalWork)
    std::cout << "*";
  std::cout << "F_tol is doing negative work against external force field by "
            << previousEnergy.externalWork - totalForceEnergy.externalWork
            << std::endl;

  // ==========================================================
  // ================     Kinetic Energy     ==================
  // ==========================================================
  std::cout << "\n";
  if (totalForceEnergy.kineticEnergy > previousEnergy.kineticEnergy)
    std::cout << "*";
  std::cout << "With F_tol, kE has increased "
            << totalForceEnergy.kineticEnergy - previousEnergy.kineticEnergy
            << " from " << previousEnergy.kineticEnergy << " to "
            << totalForceEnergy.kineticEnergy << std::endl;
}

void System::saveRichData(std::string PathToSave, bool isJustGeometry) {

  if (isJustGeometry) {
    gcs::writeSurfaceMesh(*mesh, *vpg, PathToSave);
  } else {
    gcs::RichSurfaceMeshData richData(*mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*vpg);

    // write protein distribution
    richData.addVertexProperty("protein_density", proteinDensity);
    richData.addVertexProperty("velocity", forces.ontoNormal(velocity));

    // write bool
    gcs::VertexData<double> msk(*mesh);
    msk.fromVector(toMatrix(forces.forceMask).rowwise().sum());
    richData.addVertexProperty("force_mask", msk);
    richData.addVertexProperty("protein_mask", forces.proteinMask);
    gcs::VertexData<int> tkr(*mesh);
    tkr.fromVector(centerTracker.raw().cast<int>());
    richData.addVertexProperty("the_point", tkr);
    // gcs::VertexData<int> mutMkr(*mesh);
    // mutMkr.fromVector(mutationMarker.raw().cast<int>());
    // richData.addVertexProperty("smoothing_mask", mutMkr);

    // write geometry
    gcs::VertexData<double> meanCurv(*mesh);
    meanCurv.fromVector(vpg->vertexMeanCurvatures.raw().array() /
                        vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("mean_curvature", meanCurv);
    gcs::VertexData<double> gaussCurv(*mesh);
    gaussCurv.fromVector(vpg->vertexGaussianCurvatures.raw().array() /
                         vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("gauss_curvature", gaussCurv);
    richData.addVertexProperty("spon_curvature", H0);
    richData.addVertexProperty("dual_area", vpg->vertexDualAreas);

    // write pressures
    richData.addVertexProperty("bending_force", forces.bendingForce);
    richData.addVertexProperty("deviatoric_force", forces.deviatoricForce);
    richData.addVertexProperty("capillary_force", forces.capillaryForce);
    richData.addVertexProperty("line_tension_force", forces.lineCapillaryForce);
    richData.addVertexProperty("osmotic_force", forces.osmoticForce);
    richData.addVertexProperty("adsorption_force", forces.adsorptionForce);
    richData.addVertexProperty("aggregation_force", forces.aggregationForce);
    richData.addVertexProperty("external_force", forces.externalForce);
    richData.addVertexProperty("avoidance_force", forces.selfAvoidanceForce);
    richData.addVertexProperty("physical_force", forces.mechanicalForce);

    // write chemical potential
    richData.addVertexProperty("diffusion_potential",
                               forces.diffusionPotential);
    richData.addVertexProperty("bending_potential", forces.bendingPotential);
    richData.addVertexProperty("deviatoric_potential",
                               forces.deviatoricPotential);
    richData.addVertexProperty("adsorption_potential",
                               forces.adsorptionPotential);
    richData.addVertexProperty("aggregation_potential",
                               forces.aggregationPotential);
    richData.addVertexProperty("chemical_potential", forces.chemicalPotential);

    richData.write(PathToSave);
  }
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
      mem3dg_runtime_message("Velocity is not finite!");
    }

    if (!std::isfinite(toMatrix(forces.mechanicalForceVec).norm())) {
      if (!std::isfinite(toMatrix(forces.capillaryForceVec).norm())) {
        mem3dg_runtime_message("Capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.adsorptionForceVec).norm())) {
        mem3dg_runtime_message("Adsorption force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.aggregationForceVec).norm())) {
        mem3dg_runtime_message("Aggregation force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.bendingForceVec).norm())) {
        mem3dg_runtime_message("Bending force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.deviatoricForceVec).norm())) {
        mem3dg_runtime_message("Deviatoric force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.osmoticForceVec).norm())) {
        mem3dg_runtime_message("Osmotic force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.lineCapillaryForceVec).norm())) {
        mem3dg_runtime_message("Line capillary force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.externalForceVec).norm())) {
        mem3dg_runtime_message("External force is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.selfAvoidanceForceVec).norm())) {
        mem3dg_runtime_message("Self avoidance force is not finite!");
      }
    }
  }

  if (!std::isfinite(chemErrorNorm)) {
    finite = false;

    if (!std::isfinite(toMatrix(proteinVelocity).norm())) {
      mem3dg_runtime_message("Protein velocity is not finite!");
    }

    if (!std::isfinite(toMatrix(forces.chemicalPotential).norm())) {
      if (!std::isfinite(toMatrix(forces.bendingPotential).norm())) {
        mem3dg_runtime_message("Bending Potential is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.deviatoricPotential).norm())) {
        mem3dg_runtime_message("Deviatoric Potential is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.interiorPenaltyPotential).norm())) {
        mem3dg_runtime_message(
            "Protein interior penalty potential is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.diffusionPotential).norm())) {
        mem3dg_runtime_message("Diffusion potential is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.adsorptionPotential).norm())) {
        mem3dg_runtime_message("Adsorption potential is not finite!");
      }
      if (!std::isfinite(toMatrix(forces.aggregationPotential).norm())) {
        mem3dg_runtime_message("Aggregation potential is not finite!");
      }
    }
  }

  if (!std::isfinite(energy.totalEnergy)) {
    finite = false;
    if (!std::isfinite(energy.kineticEnergy)) {
      mem3dg_runtime_message("Kinetic energy is not finite!");
    }
    if (!std::isfinite(energy.externalWork)) {
      mem3dg_runtime_message("External work is not finite!");
    }
    if (!std::isfinite(energy.potentialEnergy)) {
      if (!std::isfinite(energy.bendingEnergy)) {
        mem3dg_runtime_message("Bending energy is not finite!");
      }
      if (!std::isfinite(energy.deviatoricEnergy)) {
        mem3dg_runtime_message("Deviatoric energy is not finite!");
      }
      if (!std::isfinite(energy.surfaceEnergy)) {
        mem3dg_runtime_message("Surface energy is not finite!");
      }
      if (!std::isfinite(energy.pressureEnergy)) {
        mem3dg_runtime_message("Pressure energy is not finite!");
      }
      if (!std::isfinite(energy.adsorptionEnergy)) {
        mem3dg_runtime_message("Adsorption energy is not finite!");
      }
      if (!std::isfinite(energy.aggregationEnergy)) {
        mem3dg_runtime_message("Aggregation energy is not finite!");
      }
      if (!std::isfinite(energy.dirichletEnergy)) {
        mem3dg_runtime_message("Line tension energy is not finite!");
      }
      if (!std::isfinite(energy.proteinInteriorPenalty)) {
        mem3dg_runtime_message(
            "Protein interior penalty energy is not finite!");
      }
      if (!std::isfinite(energy.selfAvoidancePenalty)) {
        mem3dg_runtime_message(
            "Membrane self-avoidance penalty energy is not finite!");
      }
    }
  }

  return finite;
}

double System::inferTargetSurfaceArea() {
  double targetArea;
  if (isOpenMesh) {
    targetArea = parameters.tension.A_res;
    for (gcs::BoundaryLoop bl : mesh->boundaryLoops()) {
      targetArea += computePolygonArea(bl, vpg->inputVertexPositions);
    }
  } else {
    targetArea = vpg->faceAreas.raw().sum();
  }
  return targetArea;
}

void System::findVertexCenter(gcs::VertexPositionGeometry &vpg,
                              gcs::VertexData<double> &geodesicDistance,
                              double range) {
  if (parameters.point.isFloatVertex)
    mem3dg_runtime_error("parameters.point.isFloatVertex is activated!");

  if (parameters.point.pt.rows() == 1) {
    // Assign surface point as the indexed vertex
    center =
        gc::SurfacePoint(mesh->vertex((std::size_t)parameters.point.pt[0]));

  } else if (parameters.point.pt.rows() == 2 ||
             parameters.point.pt.rows() == 3) {
    // Find the cloest vertex to the point in the x-y plane or in the space
    center = gc::SurfacePoint(closestVertexToPt(*mesh, vpg, parameters.point.pt,
                                                geodesicDistance, range));
  } else {
    mem3dg_runtime_error("parameters.point.pt type not supported!");
  }

  centerTracker.fill(false);
  centerTracker[center.vertex] = true;
  int num = centerTracker.raw().cast<int>().sum();
  if (num == 0) {
    mem3dg_runtime_error("Surface point is not updated!");
  } else if (num > 1) {
    mem3dg_runtime_error("there is no "
                         "unique/existing center!");
  }
}

void System::findFloatCenter(gcs::VertexPositionGeometry &vpg,
                             gcs::VertexData<double> &geodesicDistance,
                             double range) {
  if (!parameters.point.isFloatVertex)
    mem3dg_runtime_error("parameters.point.isFloatVertex is not activated!");

  if (parameters.point.pt.rows() == 1) {
    throw std::logic_error(
        "To have Floating vertex, one must specify vertex by coordinate!");
  } else if (parameters.point.pt.rows() == 2) {
    // Find the cloest vertex to the point in the x-y plane
    gcs::Vertex closestVertex = closestVertexToPt(
        *mesh, vpg, parameters.point.pt, geodesicDistance, range);
    double shortestDistance = 1e18;
    // loop over every faces around the vertex
    for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
      if (he.isInterior()) {
        // specify vertex coordinates and the target coordinate on the face
        gc::Vector2 v1{vpg.inputVertexPositions[he.vertex()].x,
                       vpg.inputVertexPositions[he.vertex()].y};
        gc::Vector2 v2{vpg.inputVertexPositions[he.next().vertex()].x,
                       vpg.inputVertexPositions[he.next().vertex()].y};
        gc::Vector2 v3{vpg.inputVertexPositions[he.next().next().vertex()].x,
                       vpg.inputVertexPositions[he.next().next().vertex()].y};
        gc::Vector2 v{parameters.point.pt[0], parameters.point.pt[1]};
        // find the inverse barycentric mapping based on the cartesian
        // vertex coordinates
        gc::Vector3 baryCoords_ = cartesianToBarycentric(v1, v2, v3, v);

        if (baryCoords_.x > 0 && baryCoords_.y > 0 &&
            baryCoords_.z > 0) { // A. set the surface point when the point
                                 // lays within the triangle
          center = gcs::SurfacePoint(
              he.face(), correspondBarycentricCoordinates(baryCoords_, he));
          break;
        } else { // B. avoid the floating point comparision, find the best
                 // by looping over the whole fan
          baryCoords_ = gc::componentwiseMax(baryCoords_, gc::Vector3{0, 0, 0});
          baryCoords_ /= gc::sum(baryCoords_);
          gcs::SurfacePoint someSurfacePoint(
              he.face(), correspondBarycentricCoordinates(baryCoords_, he));
          double distance =
              (gc::Vector2{parameters.point.pt[0], parameters.point.pt[1]} -
               gc::Vector2{
                   someSurfacePoint.interpolate(vpg.inputVertexPositions).x,
                   someSurfacePoint.interpolate(vpg.inputVertexPositions).y})
                  .norm();
          if (distance < shortestDistance) {
            center = someSurfacePoint;
            shortestDistance = distance;
          }
        }
      }
    }
  } else if (parameters.point.pt.rows() == 3) {
    // initialize embedded point and the closest vertex
    gc::Vector3 embeddedPoint{parameters.point.pt[0], parameters.point.pt[1],
                              parameters.point.pt[2]};
    gcs::Vertex closestVertex = closestVertexToPt(
        *mesh, vpg, parameters.point.pt, geodesicDistance, range);
    gc::Vector3 vertexToPoint =
        embeddedPoint - vpg.inputVertexPositions[closestVertex];
    // initialize the surface point as the closest vertex
    // center = gc::SurfacePoint(closestVertex);
    // double shortestDistance = vertexToPoint.norm();
    double shortestDistance = 1e18;
    // loop over every faces around the vertex
    for (gcs::Halfedge he : closestVertex.outgoingHalfedges()) {
      if (he.isInterior()) {
        // project the embedded point onto the face
        auto faceNormal = vpg.faceNormal(he.face());
        gc::Vector3 projectedEmbeddedPoint =
            embeddedPoint - gc::dot(vertexToPoint, faceNormal) * faceNormal;
        // determine the choice of coordinates used for inverse
        // barycentric mapping based on orientation of the face
        gc::Vector2 v1, v2, v3, v;
        if (abs(faceNormal.z) > std::sqrt(3) / 3) {
          v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].x,
                           vpg.inputVertexPositions[he.vertex()].y};
          v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].x,
                           vpg.inputVertexPositions[he.next().vertex()].y};
          v3 = gc::Vector2{
              vpg.inputVertexPositions[he.next().next().vertex()].x,
              vpg.inputVertexPositions[he.next().next().vertex()].y};
          v = gc::Vector2{projectedEmbeddedPoint.x, projectedEmbeddedPoint.y};
        } else if (abs(faceNormal.x) > std::sqrt(3) / 3) {
          v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].y,
                           vpg.inputVertexPositions[he.vertex()].z};
          v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].y,
                           vpg.inputVertexPositions[he.next().vertex()].z};
          v3 = gc::Vector2{
              vpg.inputVertexPositions[he.next().next().vertex()].y,
              vpg.inputVertexPositions[he.next().next().vertex()].z};
          v = gc::Vector2{projectedEmbeddedPoint.y, projectedEmbeddedPoint.z};
        } else {
          v1 = gc::Vector2{vpg.inputVertexPositions[he.vertex()].z,
                           vpg.inputVertexPositions[he.vertex()].x};
          v2 = gc::Vector2{vpg.inputVertexPositions[he.next().vertex()].z,
                           vpg.inputVertexPositions[he.next().vertex()].x};
          v3 = gc::Vector2{
              vpg.inputVertexPositions[he.next().next().vertex()].z,
              vpg.inputVertexPositions[he.next().next().vertex()].x};
          v = gc::Vector2{projectedEmbeddedPoint.z, projectedEmbeddedPoint.x};
        }
        // find the inverse barycentric mapping based on the cartesian
        // vertex coordinates
        gc::Vector3 baryCoords_ = cartesianToBarycentric(v1, v2, v3, v);
        // since might not find the perfect reflecting face, best we could
        // do within each triangle
        baryCoords_ = gc::componentwiseMax(baryCoords_, gc::Vector3{0, 0, 0});
        baryCoords_ /= gc::sum(baryCoords_);
        gcs::SurfacePoint someSurfacePoint(
            he.face(), correspondBarycentricCoordinates(baryCoords_, he));
        // compute optimum distance and set surface point
        double distance = (embeddedPoint - someSurfacePoint.interpolate(
                                               vpg.inputVertexPositions))
                              .norm();
        if (distance < shortestDistance) {
          center = someSurfacePoint;
          shortestDistance = distance;
        }
      }
    }
  }

  // mark three vertex on the face
  centerTracker.fill(false);
  centerTracker[center.face.halfedge().vertex()] = true;
  centerTracker[center.face.halfedge().next().vertex()] = true;
  centerTracker[center.face.halfedge().next().next().vertex()] = true;
  int num = centerTracker.raw().cast<int>().sum();
  if (num == 0) {
    mem3dg_runtime_error("Surface point is not updated!");
  } else if (num > 3) {
    mem3dg_runtime_error("there is no "
                         "unique/existing center!");
  }
}

} // namespace solver
} // namespace mem3dg