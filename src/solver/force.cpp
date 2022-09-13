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

void System::computeGeometricForces() {
  assert(mesh->isCompressed());
  // if(!mesh->isCompressed()){
  //   mem3dg_runtime_error("Mesh must be compressed to compute forces!");
  // }

  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
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
  gc::Vertex v{mesh->vertex(i)};
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
  double Hi = vpg->vertexMeanCurvatures[i] / vpg->vertexDualAreas[i];
  double KGi = vpg->vertexGaussianCurvatures[i];
  double H0i = H0[i];
  double Kbi = Kb[i];
  double Kdi = Kd[i];
  double proteinDensityi = proteinDensity[i];
  double areaDifferenceK =
      0.25 * parameters.bending.alpha * parameters.bending.Kb * constants::PI;
  double totalMeanCurvature = vpg->vertexMeanCurvatures.raw().sum();

  bool boundaryVertex = v.isBoundary();

  for (gc::Halfedge he : v.outgoingHalfedges()) {
    std::size_t fID = he.face().getIndex();

    // Initialize local variables for computation
    std::size_t i_vj = he.tipVertex().getIndex();

    gc::Vector3 dphi_ijk{he.isInterior() ? proteinDensityGradient[fID]
                                         : gc::Vector3{0, 0, 0}};
    double Hj = vpg->vertexMeanCurvatures[i_vj] / vpg->vertexDualAreas[i_vj];
    double KGj = vpg->vertexGaussianCurvatures[i_vj];
    double H0j = H0[i_vj];
    double Kbj = Kb[i_vj];
    double Kdj = Kd[i_vj];
    double proteinDensityj = proteinDensity[i_vj];
    bool interiorHalfedge = he.isInterior();
    bool boundaryEdge = he.edge().isBoundary();
    bool boundaryNeighborVertex = he.next().vertex().isBoundary();

    // compute fundamental variational vectors
    gc::Vector3 areaGrad = 2 * computeHalfedgeMeanCurvatureVector(*vpg, he);
    gc::Vector3 gaussVec = computeHalfedgeGaussianCurvatureVector(*vpg, he);
    gc::Vector3 schlafliVec1 = vpg->edgeLengths[he.edge()] *
                               computeDihedralAngleVariation(he, he.vertex());
    gc::Vector3 schlafliVec2 =
        vpg->edgeLengths[he.twin().edge()] *
            computeDihedralAngleVariation(he.twin(), he.vertex()) +
        vpg->edgeLengths[he.next().edge()] *
            computeDihedralAngleVariation(he.next(), he.vertex()) +
        vpg->edgeLengths[he.twin().next().next().edge()] *
            computeDihedralAngleVariation(he.twin().next().next(), he.vertex());
    gc::Vector3 totalMeanCurvatureGrad =
        gaussVec + 0.5 * schlafliVec1 + 0.5 * schlafliVec2;
    gc::Vector3 oneSidedAreaGrad =
        interiorHalfedge ? 0.5 * gc::cross(vpg->faceNormals[fID],
                                           vecFromHalfedge(he.next(), *vpg))
                         : gc::Vector3{0.0, 0.0, 0.0};
    gc::Vector3 dirichletVec =
        interiorHalfedge
            ? computeHalfedgeSquaredIntegratedDerivativeNormVariationVector(
                  proteinDensity, he) /
                  vpg->faceAreas[fID]
            : gc::Vector3{0.0, 0.0, 0.0};

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
      osmoticForceVec += forces.osmoticPressure *
                         computeHalfedgeVolumeVariationVector(*vpg, he);
    if (forces.surfaceTension != 0) // surface capillary force
      capillaryForceVec -= forces.surfaceTension * areaGrad;
    if (Kdi != 0 || Kdj != 0) { // deviatoric curvature force
      if (boundaryVertex) {
        if (!boundaryEdge)
          deviatoricCurvatureForceVec_gauss +=
              2 * Kdj * KGj *
                  computeCornerAngleVariation(he.next().corner(), he.vertex()) +
              2 * Kdj * KGj *
                  computeCornerAngleVariation(he.twin().corner(), he.vertex());
      } else {
        if (boundaryNeighborVertex) {
          deviatoricCurvatureForceVec_gauss +=
              2 * Kdi * KGi *
              computeCornerAngleVariation(he.corner(), he.vertex());
        } else {
          deviatoricCurvatureForceVec_gauss +=
              2 * Kdi * KGi *
                  computeCornerAngleVariation(he.corner(), he.vertex()) +
              2 * Kdj * KGj *
                  computeCornerAngleVariation(he.next().corner(), he.vertex()) +
              2 * Kdj * KGj *
                  computeCornerAngleVariation(he.twin().corner(), he.vertex());
        }
      }
      // deviatoricCurvatureForceVec_mean -=
      //     (Kdi * Hi + Kdj * Hj) * gaussVec +
      //     (Kdi * (-Hi * Hi) / 3 + Kdj * (-Hj * Hj) * 2 / 3) * areaGrad +
      //     (Kdi * Hi * schlafliVec1 + Kdj * Hj * schlafliVec2);
      // if (boundaryVertex) {
      //   if (!boundaryEdge)
      //     deviatoricCurvatureForceVec_gauss -=
      //         Kdj *
      //             computeCornerAngleVariation(he.next().corner(),
      //             he.vertex()) +
      //         Kdj *
      //             computeCornerAngleVariation(he.twin().corner(),
      //             he.vertex());
      // } else {
      //   if (boundaryNeighborVertex) {
      //     deviatoricCurvatureForceVec_gauss -=
      //         Kdi * computeCornerAngleVariation(he.corner(), he.vertex());
      //   } else {
      //     deviatoricCurvatureForceVec_gauss -=
      //         Kdi * computeCornerAngleVariation(he.corner(), he.vertex()) +
      //         Kdj *
      //             computeCornerAngleVariation(he.next().corner(),
      //             he.vertex()) +
      //         Kdj *
      //             computeCornerAngleVariation(he.twin().corner(),
      //             he.vertex());
      //   }
      // }
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
              (2 * totalMeanCurvature * totalMeanCurvatureGrad * surfaceArea -
               totalMeanCurvature * totalMeanCurvature * areaGrad) /
              surfaceArea / surfaceArea -
          4 * parameters.bending.dA0 * areaDifferenceK / parameters.bending.D *
              (totalMeanCurvatureGrad * surfaceArea -
               totalMeanCurvature * areaGrad) /
              surfaceArea / surfaceArea -
          parameters.bending.dA0 * parameters.bending.dA0 * areaDifferenceK /
              parameters.bending.D / parameters.bending.D / surfaceArea /
              surfaceArea * areaGrad;
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
    toMatrix(forces.externalForceVec) =
        forces.maskForce(parameters.external.form(
            toMatrix(vpg->inputVertexPositions), vpg->vertexDualAreas.raw(),
            time, geodesicDistance.raw()));
    forces.externalForce = forces.ontoNormal(forces.externalForceVec);
  }
  return toMatrix(forces.externalForceVec);
}

void System::computeSelfAvoidanceForce() {
  forces.selfAvoidanceForceVec.fill({0, 0, 0});
  const double d0 = parameters.selfAvoidance.d;
  const double mu = parameters.selfAvoidance.mu;
  const double n = parameters.selfAvoidance.n;
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex vi{mesh->vertex(i)};
    gc::VertexData<bool> neighborList(*mesh, false);
    meshProcessor.meshMutator.markVertices(neighborList, vi, n);
    for (std::size_t j = i + 1; j < mesh->nVertices(); ++j) {
      if (neighborList[j])
        continue;
      gc::Vertex vj{mesh->vertex(j)};
      // double penalty = mu * vpg->vertexDualAreas[vi] * proteinDensity[vi] *
      //                  vpg->vertexDualAreas[vj] * proteinDensity[vj];
      double penalty = mu * proteinDensity[vi] * proteinDensity[vj];
      // double penalty = mu;
      // double penalty = mu * vpg->vertexDualAreas[vi] *
      // vpg->vertexDualAreas[vj];;
      gc::Vector3 r =
          vpg->inputVertexPositions[vj] - vpg->inputVertexPositions[vi];
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
  gcs::VertexData<double> dH0dphi(*mesh, 0);
  gcs::VertexData<double> dKbdphi(*mesh, 0);
  gcs::VertexData<double> dKddphi(*mesh, 0);
  auto meanCurvDiff = (vpg->vertexMeanCurvatures.raw().array() /
                       vpg->vertexDualAreas.raw().array()) -
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
        -vpg->vertexDualAreas.raw().array() *
        (meanCurvDiff * meanCurvDiff * dKbdphi.raw().array() -
         2 * Kb.raw().array() * meanCurvDiff * dH0dphi.raw().array()));
  }

  if (parameters.bending.Kd != 0 || parameters.bending.Kdc != 0) {
    forces.deviatoricCurvaturePotential.raw() =
        -dKddphi.raw().array() *
        vpg->vertexGaussianCurvatures.raw().array().square();
    // (vpg->vertexMeanCurvatures.raw().array().square() /
    //      vpg->vertexDualAreas.raw().array() -
    //  vpg->vertexGaussianCurvatures.raw().array());
  }

  if (parameters.adsorption.epsilon != 0)
    forces.adsorptionPotential.raw() = forces.maskProtein(
        -parameters.adsorption.epsilon * vpg->vertexDualAreas.raw().array());

  if (parameters.aggregation.chi != 0)
    forces.aggregationPotential.raw() = forces.maskProtein(
        -32 * parameters.aggregation.chi * (proteinDensity.raw().array() - 1) *
        (2 * proteinDensity.raw().array() - 1) * proteinDensity.raw().array() *
        vpg->vertexDualAreas.raw().array());

  if (parameters.entropy.xi != 0) {
    forces.entropyPotential.raw() =
        forces.maskProtein(-parameters.entropy.xi *
                           (proteinDensity.raw().array().log() -
                            (1 - proteinDensity.raw().array()).log()) *
                           vpg->vertexDualAreas.raw().array());
    // forces.entropyPotential.raw() = forces.maskProtein(
    //     -parameters.entropy.xi * proteinDensity.raw().array().log() *
    //     vpg->vertexDualAreas.raw().array());
  }

  if (parameters.dirichlet.eta != 0)
    forces.dirichletPotential.raw() = forces.maskProtein(
        -parameters.dirichlet.eta * vpg->cotanLaplacian * proteinDensity.raw());

  if (parameters.protein.proteinInteriorPenalty != 0)
    forces.interiorPenaltyPotential.raw() =
        forces.maskProtein(parameters.protein.proteinInteriorPenalty *
                           (1 / proteinDensity.raw().array() -
                            1 / (1 - proteinDensity.raw().array())));

  forces.chemicalPotential =
      forces.adsorptionPotential + forces.aggregationPotential +
      forces.entropyPotential + forces.spontaneousCurvaturePotential +
      forces.deviatoricCurvaturePotential + forces.dirichletPotential +
      forces.interiorPenaltyPotential;
}

EigenVectorX1d
System::computeInPlaneFluxForm(EigenVectorX1d &chemicalPotential) {
  gcs::EdgeData<double> edgeProteinDensity(*mesh, 0);
  for (std::size_t i = 0; i < mesh->nEdges(); ++i) {
    gc::Edge e{mesh->edge(i)};
    gc::Halfedge he = e.halfedge();
    edgeProteinDensity[i] = 0.5 * (proteinDensity[he.tailVertex()] +
                                   proteinDensity[he.tipVertex()]);
  }
  return vpg->hodge1 * edgeProteinDensity.raw().asDiagonal() * vpg->d0 *
         vpg->hodge0Inverse * chemicalPotential;
}

void System::computeDPDForces(double dt) {
  toMatrix(forces.dampingForceVec).setZero();
  toMatrix(forces.stochasticForceVec).setZero();
  // std::default_random_engine random_generator;
  // gcs::EdgeData<double> random_var(mesh);
  double sigma = sqrt(2 * parameters.dpd.gamma * mem3dg::constants::kBoltzmann *
                      parameters.temperature / dt);
  std::normal_distribution<double> normal_dist(0, sigma);

  for (gcs::Edge e : mesh->edges()) {
    gcs::Halfedge he = e.halfedge();
    gcs::Vertex v1 = he.vertex();
    gcs::Vertex v2 = he.next().vertex();

    gc::Vector3 dVel12 = velocity[v1] - velocity[v2];
    gc::Vector3 direction =
        (vpg->inputVertexPositions[v1] - vpg->inputVertexPositions[v2])
            .normalize();
    // gc::Vector3 direction =
    //     (vpg->vertexNormals[v1] + vpg->vertexNormals[v2]).normalize();

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
  for (std::size_t i = 0; i < mesh->nVertices(); ++i) {
    gc::Vertex v{mesh->vertex(i)};
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

          double lcr = computeLengthCrossRatio(*vpg, ij);
          double lcr0 = refLcrs[ij];
          gc::Vector3 grad_li = vecFromHalfedge(li, *vpg).normalize();
          gc::Vector3 grad_ik = vecFromHalfedge(ik.twin(), *vpg).normalize();
          gc::Vector3 lcrGrad = vpg->edgeLengths[kj.edge()] /
                                vpg->edgeLengths[jl.edge()] *
                                (grad_li * vpg->edgeLengths[ik.edge()] -
                                 grad_ik * vpg->edgeLengths[li.edge()]) /
                                pow(vpg->edgeLengths[ik.edge()], 2);
          lcrSpringForceVec -=
              parameters.spring.Kst * lcrGrad * (lcr - lcr0) / lcr0 / lcr0;
        }
        ij = he.next();
        if (!ij.edge().isBoundary()) {
          gcs::Halfedge jl = ij.next();
          gcs::Halfedge li = jl.next();
          gcs::Halfedge ik = ij.twin().next();
          gcs::Halfedge kj = ik.next();

          double lcr = computeLengthCrossRatio(*vpg, ij);
          double lcr0 = refLcrs[ij];
          gc::Vector3 grad_li = vecFromHalfedge(li.twin(), *vpg).normalize();
          gc::Vector3 grad_jl = vecFromHalfedge(jl, *vpg).normalize();
          gc::Vector3 lcrGrad = vpg->edgeLengths[kj.edge()] /
                                vpg->edgeLengths[ik.edge()] *
                                (grad_li * vpg->edgeLengths[jl.edge()] -
                                 grad_jl * vpg->edgeLengths[li.edge()]) /
                                pow(vpg->edgeLengths[jl.edge()], 2);
          lcrSpringForceVec -=
              parameters.spring.Kst * lcrGrad * (lcr - lcr0) / lcr0 / lcr0;
        }
      }

      // Local area regularization
      if (parameters.spring.Ksl != 0 && he.isInterior()) {
        gcs::Halfedge base_he = he.next();
        gc::Vector3 base_vec = vecFromHalfedge(base_he, *vpg);
        gc::Vector3 localAreaGradient =
            0.5 * gc::cross(vpg->faceNormal(he.face()), base_vec);
        double a0 = refVpg->faceAreas[base_he.face()];
        faceSpringForceVec -= parameters.spring.Ksl * localAreaGradient *
                              (vpg->faceAreas[base_he.face()] - a0) / a0 / a0;
      }

      // local edge regularization
      if (parameters.spring.Kse != 0) {
        double l0 = refVpg->edgeLengths[e];
        gc::Vector3 edgeGradient = -vecFromHalfedge(he, *vpg).normalize();
        edgeSpringForceVec -= parameters.spring.Kse * edgeGradient *
                              (vpg->edgeLengths[e] - l0) / l0 / l0;
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
