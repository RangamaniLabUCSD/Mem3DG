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

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#endif

#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void System::computeGeometricForces() {
  MEM3DG_SAFETY_ASSERT(geometry.mesh->isCompressed(),
                       "Mesh must be compressed to compute geometric forces.");
#pragma omp parallel for
  for (std::size_t i = 0; i < geometry.mesh->nVertices(); ++i) {
    computeGeometricForces(i);
  }

  // measure smoothness
  // if (meshProcessor.meshMutator.isSplitEdge ||
  //     meshProcessor.meshMutator.isCollapseEdge) {
  //   isSmooth = !hasOutlier(forces.spontaneousCurvatureForce.raw());
  // }
}

void System::computeGeometricForces(gcs::Vertex &v) {
  size_t i = v.getIndex();
  computeGeometricForces(i);
}

void System::computeGeometricForces(size_t i) {
  gc::Vertex v{geometry.mesh->vertex(i)};
  gc::Vector3 spontaneousCurvatureForceVec{0, 0, 0};
  gc::Vector3 spontaneousCurvatureForceVec_areaGrad{0, 0, 0};
  gc::Vector3 spontaneousCurvatureForceVec_gaussVec{0, 0, 0};
  gc::Vector3 spontaneousCurvatureForceVec_schlafliVec{0, 0, 0};
  gc::Vector3 deviatoricCurvatureForceVec{0, 0, 0};
  gc::Vector3 deviatoricCurvatureForceVec_mean{0, 0, 0};
  gc::Vector3 deviatoricCurvatureForceVec_gauss{0, 0, 0};
  gc::Vector3 areaDifferenceForceVec{0, 0, 0};
  gc::Vector3 capillaryForceVec{0, 0, 0};
  gc::Vector3 osmoticForceVec{0, 0, 0};
  gc::Vector3 lineCapillaryForceVec{0, 0, 0};
  gc::Vector3 adsorptionForceVec{0, 0, 0};
  gc::Vector3 aggregationForceVec{0, 0, 0};
  gc::Vector3 entropyForceVec{0, 0, 0};
  double Ai = geometry.vpg->vertexDualAreas[i];
  double Hi = geometry.vpg->vertexMeanCurvatures[i] / Ai;
  double KGi = geometry.vpg->vertexGaussianCurvatures[i] / Ai;
  double H0i = H0[i];
  double Kbi = Kb[i];
  double Kdi = Kd[i];
  double proteinDensityi = proteinDensity[i];
  double areaDifferenceK =
      0.25 * parameters.bending.alpha * parameters.bending.Kb * constants::PI;
  double totalMeanCurvature = geometry.vpg->vertexMeanCurvatures.raw().sum();

  bool boundaryVertex = v.isBoundary();

  for (gc::Halfedge he : v.outgoingHalfedges()) {
    std::size_t fID = he.face().getIndex();

    // Initialize local variables for computation
    std::size_t i_vj = he.tipVertex().getIndex();

    gc::Vector3 dphi_ijk{he.isInterior() ? proteinDensityGradient[fID]
                                         : gc::Vector3{0, 0, 0}};
    double Aj = geometry.vpg->vertexDualAreas[i_vj];
    double Hj = geometry.vpg->vertexMeanCurvatures[i_vj] / Aj;
    double KGj = geometry.vpg->vertexGaussianCurvatures[i_vj] / Aj;
    double H0j = H0[i_vj];
    double Kbj = Kb[i_vj];
    double Kdj = Kd[i_vj];
    double proteinDensityj = proteinDensity[i_vj];
    bool interiorHalfedge = he.isInterior();
    bool boundaryEdge = he.edge().isBoundary();
    bool boundaryNeighborVertex = he.next().vertex().isBoundary();

    // compute fundamental variational vectors
    gc::Vector3 areaGrad =
        2 * geometry.computeHalfedgeMeanCurvatureVector(*geometry.vpg, he);
    gc::Vector3 gaussVec =
        geometry.computeHalfedgeGaussianCurvatureVector(*geometry.vpg, he);
    gc::Vector3 schlafliVec1 =
        geometry.vpg->edgeLengths[he.edge()] *
        geometry.computeDihedralAngleVariation(he, he.vertex());
    gc::Vector3 schlafliVec2 =
        geometry.vpg->edgeLengths[he.twin().edge()] *
            geometry.computeDihedralAngleVariation(he.twin(), he.vertex()) +
        geometry.vpg->edgeLengths[he.next().edge()] *
            geometry.computeDihedralAngleVariation(he.next(), he.vertex()) +
        geometry.vpg->edgeLengths[he.twin().next().next().edge()] *
            geometry.computeDihedralAngleVariation(he.twin().next().next(),
                                                   he.vertex());
    gc::Vector3 totalMeanCurvatureGrad =
        gaussVec + 0.5 * schlafliVec1 + 0.5 * schlafliVec2;
    gc::Vector3 oneSidedAreaGrad =
        interiorHalfedge
            ? 0.5 * gc::cross(geometry.vpg->faceNormals[fID],
                              vecFromHalfedge(he.next(), *geometry.vpg))
            : gc::Vector3{0.0, 0.0, 0.0};
    gc::Vector3 dirichletVec =
        interiorHalfedge
            ? geometry.computeHalfedgeSquaredIntegratedDerivativeNormVariationVector(
                  proteinDensity, he) /
                  geometry.vpg->faceAreas[fID]
            : gc::Vector3{0.0, 0.0, 0.0};
    gc::Vector3 gaussVarVec1 =
        boundaryVertex ? gc::Vector3{0.0, 0.0, 0.0}
                       : -geometry.computeCornerAngleVariation(he, he.vertex());
    gc::Vector3 gaussVarVec2 =
        boundaryNeighborVertex
            ? gc::Vector3{0.0, 0.0, 0.0}
            : -geometry.computeCornerAngleVariation(he.next(), he.vertex()) -
                  geometry.computeCornerAngleVariation(he.twin(), he.vertex());

    // Assemble to forces
    if (Kbi != 0 || Kbj != 0) { // spontaneous curvature force
      spontaneousCurvatureForceVec_schlafliVec -=
          (Kbi * (Hi - H0i) * schlafliVec1 + Kbj * (Hj - H0j) * schlafliVec2);
      spontaneousCurvatureForceVec_areaGrad -=
          (Kbi * (H0i * H0i - Hi * Hi) / 3 +
           Kbj * (H0j * H0j - Hj * Hj) * 2 / 3) *
          areaGrad;
      spontaneousCurvatureForceVec_gaussVec -=
          (Kbi * (Hi - H0i) + Kbj * (Hj - H0j)) * gaussVec;
    }
    if (forces.osmoticPressure != 0) // omostic force
      osmoticForceVec +=
          forces.osmoticPressure *
          geometry.computeHalfedgeVolumeVariationVector(*geometry.vpg, he);
    if (forces.surfaceTension != 0) // surface capillary force
      capillaryForceVec -= forces.surfaceTension * areaGrad;
    if (Kdi != 0 || Kdj != 0) { // deviatoric curvature force
      deviatoricCurvatureForceVec_gauss -=
          2 * Kdi * KGi * gaussVarVec1 + 2 * Kdj * KGj * gaussVarVec2 -
          (Kdi * KGi * KGi / 3 + Kdj * KGj * KGj * 2 / 3) * areaGrad;
      // deviatoricCurvatureForceVec_mean -=
      //     (Kdi * Hi + Kdj * Hj) * gaussVec +
      //     (Kdi * (-Hi * Hi) / 3 + Kdj * (-Hj * Hj) * 2 / 3) * areaGrad +
      //     (Kdi * Hi * schlafliVec1 + Kdj * Hj * schlafliVec2);
      // deviatoricCurvatureForceVec_gauss +=
      //     Kdi * gaussVarVec1 + Kdj * gaussVarVec2;
    }
    if (parameters.adsorption.epsilon != 0) // adsorption force
      adsorptionForceVec -= (proteinDensityi / 3 + proteinDensityj * 2 / 3) *
                            parameters.adsorption.epsilon * areaGrad;
    if (parameters.aggregation.chi != 0) // aggregation force
      aggregationForceVec -=
          (pow(pow(2 * proteinDensityi - 1, 2) - 1, 2) / 3 +
           pow(pow(2 * proteinDensityj - 1, 2) - 1, 2) * 2 / 3) *
          parameters.aggregation.chi * areaGrad;
    if (parameters.entropy.xi != 0) // entropy force
      entropyForceVec -= ((proteinDensityi * log(proteinDensityi) +
                           (1 - proteinDensityi) * log(1 - proteinDensityi)) /
                              3 +
                          (proteinDensityj * log(proteinDensityj) +
                           (1 - proteinDensityj) * log(1 - proteinDensityj)) *
                              2 / 3) *
                         parameters.entropy.xi * areaGrad;
    if (parameters.dirichlet.eta != 0) // line capillary force
      lineCapillaryForceVec -=
          parameters.dirichlet.eta *
          (0.125 * dirichletVec - 0.5 * dphi_ijk.norm2() * oneSidedAreaGrad);
    if (parameters.bending.alpha != 0) // area difference force
      areaDifferenceForceVec -=
          4 * areaDifferenceK *
              (2 * totalMeanCurvature * totalMeanCurvatureGrad *
                   geometry.surfaceArea -
               totalMeanCurvature * totalMeanCurvature * areaGrad) /
              geometry.surfaceArea / geometry.surfaceArea -
          4 * parameters.bending.dA0 * areaDifferenceK / parameters.bending.D *
              (totalMeanCurvatureGrad * geometry.surfaceArea -
               totalMeanCurvature * areaGrad) /
              geometry.surfaceArea / geometry.surfaceArea -
          parameters.bending.dA0 * parameters.bending.dA0 * areaDifferenceK /
              parameters.bending.D / parameters.bending.D /
              geometry.surfaceArea / geometry.surfaceArea * areaGrad;
  }

  spontaneousCurvatureForceVec = spontaneousCurvatureForceVec_areaGrad +
                                 spontaneousCurvatureForceVec_gaussVec +
                                 spontaneousCurvatureForceVec_schlafliVec;
  deviatoricCurvatureForceVec =
      deviatoricCurvatureForceVec_mean + deviatoricCurvatureForceVec_gauss;

  // masking
  spontaneousCurvatureForceVec_areaGrad =
      forces.maskForce(spontaneousCurvatureForceVec_areaGrad, i);
  spontaneousCurvatureForceVec_gaussVec =
      forces.maskForce(spontaneousCurvatureForceVec_gaussVec, i);
  spontaneousCurvatureForceVec_schlafliVec =
      forces.maskForce(spontaneousCurvatureForceVec_schlafliVec, i);
  spontaneousCurvatureForceVec =
      forces.maskForce(spontaneousCurvatureForceVec, i);
  deviatoricCurvatureForceVec_mean =
      forces.maskForce(deviatoricCurvatureForceVec_mean, i);
  deviatoricCurvatureForceVec_gauss =
      forces.maskForce(deviatoricCurvatureForceVec_gauss, i);
  deviatoricCurvatureForceVec =
      forces.maskForce(deviatoricCurvatureForceVec, i);
  areaDifferenceForceVec = forces.maskForce(areaDifferenceForceVec, i);
  osmoticForceVec = forces.maskForce(osmoticForceVec, i);
  capillaryForceVec = forces.maskForce(capillaryForceVec, i);
  lineCapillaryForceVec = forces.maskForce(lineCapillaryForceVec, i);
  adsorptionForceVec = forces.maskForce(adsorptionForceVec, i);
  aggregationForceVec = forces.maskForce(aggregationForceVec, i);
  entropyForceVec = forces.maskForce(entropyForceVec, i);

  // Combine to one
  forces.spontaneousCurvatureForceVec_areaGrad[i] =
      spontaneousCurvatureForceVec_areaGrad;
  forces.spontaneousCurvatureForceVec_gaussVec[i] =
      spontaneousCurvatureForceVec_gaussVec;
  forces.spontaneousCurvatureForceVec_schlafliVec[i] =
      spontaneousCurvatureForceVec_schlafliVec;
  forces.spontaneousCurvatureForceVec[i] = spontaneousCurvatureForceVec;
  forces.deviatoricCurvatureForceVec[i] = deviatoricCurvatureForceVec;
  forces.deviatoricCurvatureForceVec_mean[i] = deviatoricCurvatureForceVec_mean;
  forces.deviatoricCurvatureForceVec_gauss[i] =
      deviatoricCurvatureForceVec_gauss;
  forces.areaDifferenceForceVec[i] = areaDifferenceForceVec;
  forces.capillaryForceVec[i] = capillaryForceVec;
  forces.osmoticForceVec[i] = osmoticForceVec;
  forces.lineCapillaryForceVec[i] = lineCapillaryForceVec;
  forces.adsorptionForceVec[i] = adsorptionForceVec;
  forces.aggregationForceVec[i] = aggregationForceVec;
  forces.entropyForceVec[i] = entropyForceVec;

  // Scalar force by projection to angle-weighted normal
  forces.spontaneousCurvatureForce[i] =
      forces.ontoNormal(spontaneousCurvatureForceVec, i);
  forces.deviatoricCurvatureForce[i] =
      forces.ontoNormal(deviatoricCurvatureForceVec, i);
  forces.areaDifferenceForce[i] = forces.ontoNormal(areaDifferenceForceVec, i);
  forces.capillaryForce[i] = forces.ontoNormal(capillaryForceVec, i);
  forces.osmoticForce[i] = forces.ontoNormal(osmoticForceVec, i);
  forces.lineCapillaryForce[i] = forces.ontoNormal(lineCapillaryForceVec, i);
  forces.adsorptionForce[i] = forces.ontoNormal(adsorptionForceVec, i);
  forces.aggregationForce[i] = forces.ontoNormal(aggregationForceVec, i);
  forces.entropyForce[i] = forces.ontoNormal(entropyForceVec, i);
}

EigenVectorX3dr System::prescribeExternalForce() {
  if (parameters.external.form != NULL) {
    toMatrix(forces.externalForceVec) = forces.maskForce(
        parameters.external.form(toMatrix(geometry.vpg->inputVertexPositions),
                                 geometry.vpg->vertexDualAreas.raw(), time,
                                 geometry.geodesicDistance.raw()));
    forces.externalForce = forces.ontoNormal(forces.externalForceVec);
  }
  return toMatrix(forces.externalForceVec);
}

void System::computeSelfAvoidanceForce() {
  forces.selfAvoidanceForceVec.fill({0, 0, 0});
  const double d0 = parameters.selfAvoidance.d;
  const double mu = parameters.selfAvoidance.mu;
  const double n = parameters.selfAvoidance.n;
  for (std::size_t i = 0; i < geometry.mesh->nVertices(); ++i) {
    gc::Vertex vi{geometry.mesh->vertex(i)};
    gc::VertexData<bool> neighborList(*geometry.mesh, false);
    meshProcessor.meshMutator.markVertices(neighborList, vi, n);
    for (std::size_t j = i + 1; j < geometry.mesh->nVertices(); ++j) {
      if (neighborList[j])
        continue;
      gc::Vertex vj{geometry.mesh->vertex(j)};
      // double penalty = mu * geometry.vpg->vertexDualAreas[vi] *
      // proteinDensity[vi] *
      //                  geometry.vpg->vertexDualAreas[vj] *
      //                  proteinDensity[vj];
      double penalty = mu * proteinDensity[vi] * proteinDensity[vj];
      // double penalty = mu;
      // double penalty = mu * geometry.vpg->vertexDualAreas[vi] *
      // geometry.vpg->vertexDualAreas[vj];;
      gc::Vector3 r = geometry.vpg->inputVertexPositions[vj] -
                      geometry.vpg->inputVertexPositions[vi];
      double distance = gc::norm(r) - d0;
      gc::Vector3 grad = r.normalize();
      // forces.selfAvoidanceForceVec[i] -=
      //     forces.maskForce(penalty / distance * grad, i);
      // forces.selfAvoidanceForceVec[j] +=
      //     forces.maskForce(penalty / distance * grad, j);
      forces.selfAvoidanceForceVec[i] -=
          forces.maskForce(penalty / distance / distance * grad, i);
      forces.selfAvoidanceForceVec[j] +=
          forces.maskForce(penalty / distance / distance * grad, j);
    }
  }
  forces.selfAvoidanceForce = forces.ontoNormal(forces.selfAvoidanceForceVec);
}

void System::computeChemicalPotentials() {
  gcs::VertexData<double> dH0dphi(*geometry.mesh, 0);
  gcs::VertexData<double> dKbdphi(*geometry.mesh, 0);
  gcs::VertexData<double> dKddphi(*geometry.mesh, 0);
  auto meanCurvDiff = (geometry.vpg->vertexMeanCurvatures.raw().array() /
                       geometry.vpg->vertexDualAreas.raw().array()) -
                      H0.raw().array();

  if (parameters.bending.relation == "linear") {
    dH0dphi.fill(parameters.bending.H0c);
    dKbdphi.fill(parameters.bending.Kbc);
    dKddphi.fill(parameters.bending.Kdc);
  } else if (parameters.bending.relation == "hill") {
    EigenVectorX1d proteinDensitySq =
        (proteinDensity.raw().array() * proteinDensity.raw().array()).matrix();
    dH0dphi.raw() =
        (2 * parameters.bending.H0c * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
    dKbdphi.raw() =
        (2 * parameters.bending.Kbc * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
    dKddphi.raw() =
        (2 * parameters.bending.Kdc * proteinDensity.raw().array() /
         ((1 + proteinDensitySq.array()) * (1 + proteinDensitySq.array())))
            .matrix();
  }

  if (parameters.bending.Kb != 0 || parameters.bending.Kbc != 0) {
    forces.spontaneousCurvaturePotential.raw() = forces.maskProtein(
        -geometry.vpg->vertexDualAreas.raw().array() *
        (meanCurvDiff * meanCurvDiff * dKbdphi.raw().array() -
         2 * Kb.raw().array() * meanCurvDiff * dH0dphi.raw().array()));
  }

  if (parameters.bending.Kd != 0 || parameters.bending.Kdc != 0) {
    forces.deviatoricCurvaturePotential.raw() =
        -dKddphi.raw().array() *
        geometry.vpg->vertexGaussianCurvatures.raw().array().square() /
        geometry.vpg->vertexDualAreas.raw().array();
    // (geometry.vpg->vertexMeanCurvatures.raw().array().square() /
    //      geometry.vpg->vertexDualAreas.raw().array() -
    //  geometry.vpg->vertexGaussianCurvatures.raw().array());
  }

  if (parameters.adsorption.epsilon != 0)
    forces.adsorptionPotential.raw() =
        forces.maskProtein(-parameters.adsorption.epsilon *
                           geometry.vpg->vertexDualAreas.raw().array());

  if (parameters.aggregation.chi != 0)
    forces.aggregationPotential.raw() = forces.maskProtein(
        -32 * parameters.aggregation.chi * (proteinDensity.raw().array() - 1) *
        (2 * proteinDensity.raw().array() - 1) * proteinDensity.raw().array() *
        geometry.vpg->vertexDualAreas.raw().array());

  if (parameters.entropy.xi != 0) {
    forces.entropyPotential.raw() =
        forces.maskProtein(-parameters.entropy.xi *
                           (proteinDensity.raw().array().log() -
                            (1 - proteinDensity.raw().array()).log()) *
                           geometry.vpg->vertexDualAreas.raw().array());
    // forces.entropyPotential.raw() = forces.maskProtein(
    //     -parameters.entropy.xi * proteinDensity.raw().array().log() *
    //     geometry.vpg->vertexDualAreas.raw().array());
  }

  if (parameters.dirichlet.eta != 0)
    forces.dirichletPotential.raw() =
        forces.maskProtein(-parameters.dirichlet.eta *
                           geometry.vpg->cotanLaplacian * proteinDensity.raw());

  if (parameters.protein.proteinInteriorPenalty != 0) {
    // Doing some tricks to penalize
    EigenVectorX1d a = 1 / proteinDensity.raw().array();
    a = (a.array().isFinite()).select(a, 0);
    EigenVectorX1d b = 1 / (1 - proteinDensity.raw().array());
    b = (b.array().isFinite()).select(b, 0);
    forces.interiorPenaltyPotential.raw() = forces.maskProtein(
        parameters.protein.proteinInteriorPenalty * (a.array() - b.array()));
  }

  forces.chemicalPotential =
      forces.adsorptionPotential + forces.aggregationPotential +
      forces.entropyPotential + forces.spontaneousCurvaturePotential +
      forces.deviatoricCurvaturePotential + forces.dirichletPotential +
      forces.interiorPenaltyPotential;
}

EigenVectorX1d
System::computeInPlaneFluxForm(EigenVectorX1d &chemicalPotential) {
  gcs::EdgeData<double> edgeProteinDensity(*geometry.mesh, 0);
  for (std::size_t i = 0; i < geometry.mesh->nEdges(); ++i) {
    gc::Edge e{geometry.mesh->edge(i)};
    gc::Halfedge he = e.halfedge();
    edgeProteinDensity[i] = 0.5 * (proteinDensity[he.tailVertex()] +
                                   proteinDensity[he.tipVertex()]);
  }
  return geometry.vpg->hodge1 * edgeProteinDensity.raw().asDiagonal() *
         geometry.vpg->d0 * geometry.vpg->hodge0Inverse * chemicalPotential;
}

void System::computeDPDForces(double dt) {
  toMatrix(forces.dampingForceVec).setZero();
  toMatrix(forces.stochasticForceVec).setZero();
  // std::default_random_engine random_generator;
  // gcs::EdgeData<double> random_var(mesh);
  double sigma = sqrt(2 * parameters.dpd.gamma * mem3dg::constants::kBoltzmann *
                      parameters.temperature / dt);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : geometry.mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    gcs::Vertex v1 = he.vertex();
    gcs::Vertex v2 = he.next().vertex();

    gc::Vector3 dVel12 = velocity[v1] - velocity[v2];
    gc::Vector3 direction = (geometry.vpg->inputVertexPositions[v1] -
                             geometry.vpg->inputVertexPositions[v2])
                                .normalize();
    // gc::Vector3 direction =
    //     (geometry.vpg->vertexNormals[v1] +
    //     geometry.vpg->vertexNormals[v2]).normalize();

    gc::Vector3 df =
        parameters.dpd.gamma * (gc::dot(dVel12, direction) * direction);
    forces.dampingForceVec[v1] -= df;
    forces.dampingForceVec[v2] += df;

    if (sigma != 0) {
      double noise = normal_dist(rng);
      forces.stochasticForceVec[v1] += noise * direction;
      forces.stochasticForceVec[v2] -= noise * direction;
    }

    // gc::Vector3 dVel21 = vel[v2] - vel[v1];
    // gc::Vector3 dPos21_n = (pos[v2] - pos[v1]).normalize();

    // std::cout << -gamma * (gc::dot(dVel12, dPos12_n) * dPos12_n)
    //           << " == " << -gamma * (gc::dot(-dVel12, -dPos12_n) *
    //           -dPos12_n)
    //           << " == " << -gamma * (gc::dot(dVel21, dPos21_n) * dPos21_n)
    //           << std::endl;
  }
  forces.dampingForceVec = forces.maskForce(forces.dampingForceVec);
  forces.stochasticForceVec = forces.maskForce(forces.stochasticForceVec);
  // dampingForce_e =
  //     forces.maskForce(forces.addNormal(forces.ontoNormal(dampingForce_e)));
  // stochasticForce_e =
  //     forces.maskForce(forces.addNormal(forces.ontoNormal(stochasticForce_e)));
}

gc::VertexData<gc::Vector3> System::computeDampingForce() {
  return -parameters.damping * velocity;
}

void System::computeSpringForces() {
  for (std::size_t i = 0; i < geometry.mesh->nVertices(); ++i) {
    gc::Vertex v{geometry.mesh->vertex(i)};
    gc::Vector3 edgeSpringForceVec{0, 0, 0};
    gc::Vector3 lcrSpringForceVec{0, 0, 0};
    gc::Vector3 faceSpringForceVec{0, 0, 0};

    for (gcs::Halfedge he : v.outgoingHalfedges()) {
      gcs::Edge e = he.edge();

      // Conformal regularization
      if (parameters.spring.Kst != 0) {
        gcs::Halfedge ij = he;
        if (!ij.edge().isBoundary()) {
          gcs::Halfedge jl = ij.next();
          gcs::Halfedge li = jl.next();
          gcs::Halfedge ik = ij.twin().next();
          gcs::Halfedge kj = ik.next();

          double lcr = geometry.computeLengthCrossRatio(*geometry.vpg, ij);
          double lcr0 = geometry.refLcrs[ij];
          gc::Vector3 grad_li = vecFromHalfedge(li, *geometry.vpg).normalize();
          gc::Vector3 grad_ik =
              vecFromHalfedge(ik.twin(), *geometry.vpg).normalize();
          gc::Vector3 lcrGrad =
              geometry.vpg->edgeLengths[kj.edge()] /
              geometry.vpg->edgeLengths[jl.edge()] *
              (grad_li * geometry.vpg->edgeLengths[ik.edge()] -
               grad_ik * geometry.vpg->edgeLengths[li.edge()]) /
              pow(geometry.vpg->edgeLengths[ik.edge()], 2);
          lcrSpringForceVec -=
              parameters.spring.Kst * lcrGrad * (lcr - lcr0) / lcr0 / lcr0;
        }
        ij = he.next();
        if (!ij.edge().isBoundary()) {
          gcs::Halfedge jl = ij.next();
          gcs::Halfedge li = jl.next();
          gcs::Halfedge ik = ij.twin().next();
          gcs::Halfedge kj = ik.next();

          double lcr = geometry.computeLengthCrossRatio(*geometry.vpg, ij);
          double lcr0 = geometry.refLcrs[ij];
          gc::Vector3 grad_li =
              vecFromHalfedge(li.twin(), *geometry.vpg).normalize();
          gc::Vector3 grad_jl = vecFromHalfedge(jl, *geometry.vpg).normalize();
          gc::Vector3 lcrGrad =
              geometry.vpg->edgeLengths[kj.edge()] /
              geometry.vpg->edgeLengths[ik.edge()] *
              (grad_li * geometry.vpg->edgeLengths[jl.edge()] -
               grad_jl * geometry.vpg->edgeLengths[li.edge()]) /
              pow(geometry.vpg->edgeLengths[jl.edge()], 2);
          lcrSpringForceVec -=
              parameters.spring.Kst * lcrGrad * (lcr - lcr0) / lcr0 / lcr0;
        }
      }

      // Local area regularization
      if (parameters.spring.Ksl != 0 && he.isInterior()) {
        gcs::Halfedge base_he = he.next();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, *geometry.vpg);
        gc::Vector3 localAreaGradient =
            0.5 * gc::cross(geometry.vpg->faceNormal(he.face()), base_vec);
        double a0 = geometry.refVpg->faceAreas[base_he.face()];
        faceSpringForceVec -= parameters.spring.Ksl * localAreaGradient *
                              (geometry.vpg->faceAreas[base_he.face()] - a0) /
                              a0 / a0;
      }

      // local edge regularization
      if (parameters.spring.Kse != 0) {
        double l0 = geometry.refVpg->edgeLengths[e];
        gc::Vector3 edgeGradient =
            -vecFromHalfedge(he, *geometry.vpg).normalize();
        edgeSpringForceVec -= parameters.spring.Kse * edgeGradient *
                              (geometry.vpg->edgeLengths[e] - l0) / l0 / l0;
      }
    }

    edgeSpringForceVec = forces.maskForce(edgeSpringForceVec, i);
    lcrSpringForceVec = forces.maskForce(lcrSpringForceVec, i);
    faceSpringForceVec = forces.maskForce(faceSpringForceVec, i);

    // edgeSpringForceVec = forces.toTangent(edgeSpringForceVec, v);
    // lcrSpringForceVec = forces.toTangent(lcrSpringForceVec, v);
    // faceSpringForceVec = forces.toTangent(faceSpringForceVec, v);

    forces.edgeSpringForceVec[i] = edgeSpringForceVec;
    forces.faceSpringForceVec[i] = faceSpringForceVec;
    forces.lcrSpringForceVec[i] = lcrSpringForceVec;
    forces.springForceVec[i] =
        edgeSpringForceVec + faceSpringForceVec + lcrSpringForceVec;
  }
}

void System::addNonconservativeForcing(double timeStep) {

  // compute nonconservative forcing and append it to mechanicalForce(Vec)
  forces.externalForceVec.fill({0, 0, 0});
  forces.dampingForceVec.fill({0, 0, 0});
  forces.stochasticForceVec.fill({0, 0, 0});
  if (parameters.variation.isShapeVariation) {
    // compute (potentially) time dependent forces
    if (parameters.external.form != NULL) {
      prescribeExternalForce();
      forces.mechanicalForceVec += forces.externalForceVec;
    }
    // compute dissipative force that depends on velocity
    if (parameters.damping != 0)
      forces.mechanicalForceVec += computeDampingForce();
    // compute DPD force that depends on velocity and time step
    if (parameters.variation.isShapeVariation && parameters.dpd.gamma != 0) {
      computeDPDForces(timeStep);
      forces.mechanicalForceVec +=
          forces.dampingForceVec + forces.stochasticForceVec;
    }
  }
  forces.mechanicalForce = forces.ontoNormal(forces.mechanicalForceVec);
}

void System::computeConservativeForcing() {
  // zero all forcing
  forces.mechanicalForceVec.fill({0, 0, 0});
  forces.conservativeForceVec.fill({0, 0, 0});

  forces.spontaneousCurvatureForceVec.fill({0, 0, 0});
  forces.spontaneousCurvatureForceVec_areaGrad.fill({0, 0, 0});
  forces.spontaneousCurvatureForceVec_gaussVec.fill({0, 0, 0});
  forces.spontaneousCurvatureForceVec_schlafliVec.fill({0, 0, 0});

  forces.deviatoricCurvatureForceVec.fill({0, 0, 0});
  forces.areaDifferenceForceVec.fill({0, 0, 0});

  forces.capillaryForceVec.fill({0, 0, 0});
  forces.osmoticForceVec.fill({0, 0, 0});
  forces.lineCapillaryForceVec.fill({0, 0, 0});
  forces.adsorptionForceVec.fill({0, 0, 0});
  forces.aggregationForceVec.fill({0, 0, 0});
  forces.entropyForceVec.fill({0, 0, 0});
  forces.selfAvoidanceForceVec.fill({0, 0, 0});

  forces.springForceVec.fill({0, 0, 0});
  forces.faceSpringForceVec.fill({0, 0, 0});
  forces.edgeSpringForceVec.fill({0, 0, 0});
  forces.lcrSpringForceVec.fill({0, 0, 0});

  forces.mechanicalForce.raw().setZero();
  forces.conservativeForce.raw().setZero();
  forces.spontaneousCurvatureForce.raw().setZero();
  forces.deviatoricCurvatureForce.raw().setZero();
  forces.areaDifferenceForce.raw().setZero();
  forces.capillaryForce.raw().setZero();
  forces.lineCapillaryForce.raw().setZero();
  forces.externalForce.raw().setZero();
  forces.adsorptionForce.raw().setZero();
  forces.aggregationForce.raw().setZero();
  forces.entropyForce.raw().setZero();
  forces.osmoticForce.raw().setZero();
  forces.selfAvoidanceForce.raw().setZero();

  forces.chemicalPotential.raw().setZero();

  forces.dirichletPotential.raw().setZero();
  forces.spontaneousCurvaturePotential.raw().setZero();
  forces.deviatoricCurvaturePotential.raw().setZero();
  forces.adsorptionPotential.raw().setZero();
  forces.aggregationPotential.raw().setZero();
  forces.entropyPotential.raw().setZero();
  forces.interiorPenaltyPotential.raw().setZero();

  if (parameters.variation.isShapeVariation) {
    computeGeometricForces();

    // compute conservative nongeometric forces
    if (parameters.selfAvoidance.mu != 0) {
      computeSelfAvoidanceForce();
    }
    if (parameters.spring.Kse != 0 || parameters.spring.Ksl != 0 ||
        parameters.spring.Kst != 0) {
      computeSpringForces();
    }

    // summerize all conservative forces
    forces.conservativeForceVec =
        forces.osmoticForceVec + forces.capillaryForceVec +
        forces.spontaneousCurvatureForceVec +
        forces.deviatoricCurvatureForceVec + forces.areaDifferenceForceVec +
        forces.lineCapillaryForceVec + forces.adsorptionForceVec +
        forces.aggregationForceVec + forces.entropyForceVec +
        forces.selfAvoidanceForceVec + forces.springForceVec;
    forces.conservativeForce = forces.ontoNormal(forces.conservativeForceVec);

    // mechanical force includes all conservative forces
    forces.mechanicalForceVec = forces.conservativeForceVec;
    forces.mechanicalForce = forces.ontoNormal(forces.mechanicalForceVec);
  }

  // total chemical potential is summed inside function call
  if (parameters.variation.isProteinVariation) {
    computeChemicalPotentials();
  }
}

} // namespace solver
} // namespace mem3dg
