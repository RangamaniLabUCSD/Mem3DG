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

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {
namespace solver {

void System::saveRichData(std::string PathToSave, bool isJustGeometry) {

  if (isJustGeometry) {
    gcs::writeSurfaceMesh(*mesh, *vpg, PathToSave);
  } else {
    gcs::RichSurfaceMeshData richData(*mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*vpg);

    // write protein distribution
    richData.addVertexProperty("proteinDensity", proteinDensity);
    richData.addVertexProperty("velocity", forces.ontoNormal(velocity));

    // write bool
    gcs::VertexData<double> msk(*mesh);
    msk.fromVector(toMatrix(forces.forceMask).rowwise().sum());
    richData.addVertexProperty("forceMask", msk);
    richData.addVertexProperty("proteinMask", forces.proteinMask);
    gcs::VertexData<double> tkr(*mesh);
    tkr.fromVector(notableVertex.raw().cast<double>());
    richData.addVertexProperty("notableVertex", tkr);
    // gcs::VertexData<int> mutMkr(*mesh);
    // mutMkr.fromVector(mutationMarker.raw().cast<int>());
    // richData.addVertexProperty("smoothing_mask", mutMkr);

    // write geometry
    gcs::VertexData<double> meanCurv(*mesh);
    meanCurv.fromVector(vpg->vertexMeanCurvatures.raw().array() /
                        vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("meanCurvature", meanCurv);
    gcs::VertexData<double> gaussCurv(*mesh);
    gaussCurv.fromVector(vpg->vertexGaussianCurvatures.raw().array() /
                         vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("gaussianCurvature", gaussCurv);
    richData.addVertexProperty("spontaneousCurvature", H0);
    richData.addVertexProperty("dualArea", vpg->vertexDualAreas);

    // write pressures
    richData.addVertexProperty("spontaneousCurvatureForce",
                               forces.spontaneousCurvatureForce);
    richData.addVertexProperty("deviatoricCurvatureForce",
                               forces.deviatoricCurvatureForce);
    richData.addVertexProperty("areaDifferenceForce",
                               forces.areaDifferenceForce);
    richData.addVertexProperty("capillaryForce", forces.capillaryForce);
    richData.addVertexProperty("lineCapillaryForce", forces.lineCapillaryForce);
    richData.addVertexProperty("osmoticForce", forces.osmoticForce);
    richData.addVertexProperty("adsorptionForce", forces.adsorptionForce);
    richData.addVertexProperty("aggregationForce", forces.aggregationForce);
    richData.addVertexProperty("externalForce", forces.externalForce);
    richData.addVertexProperty("selfAvoidanceForce", forces.selfAvoidanceForce);
    richData.addVertexProperty("mechanicalForce", forces.mechanicalForce);

    // write chemical potential
    richData.addVertexProperty("dirichletPotential", forces.dirichletPotential);
    richData.addVertexProperty("spontaneousCurvaturePotential",
                               forces.spontaneousCurvaturePotential);
    richData.addVertexProperty("deviatoricCurvaturePotential",
                               forces.deviatoricCurvaturePotential);
    richData.addVertexProperty("adsorptionPotential",
                               forces.adsorptionPotential);
    richData.addVertexProperty("aggregationPotential",
                               forces.aggregationPotential);
    richData.addVertexProperty("chemicalPotential", forces.chemicalPotential);

    richData.write(PathToSave);
  }
}

void System::mapContinuationVariables(std::string plyFile) {
  std::unique_ptr<gcs::SurfaceMesh> ptrMesh_local;
  std::unique_ptr<gcs::RichSurfaceMeshData> ptrRichData_local;

  // Open mesh file
  std::tie(ptrMesh_local, ptrRichData_local) =
      gcs::RichSurfaceMeshData::readMeshAndData(plyFile);

  // Check consistent topology for continuation
  if (mesh->nFaces() != ptrMesh_local->nFaces() ||
      mesh->nVertices() != ptrMesh_local->nVertices()) {
    throw std::logic_error(
        "Topology for continuation parameters mapping is not consistent!");
  } else {
    // Map continuation variables
    time = std::stod(sliceString(plyFile, "t", "_"));
    energy.time = time;
    // frame = std::stod(sliceString(plyFile, "f", "_"));
    proteinDensity =
        ptrRichData_local->getVertexProperty<double>("protein_density")
            .reinterpretTo(*mesh);
    // vel_protein =
    //     ptrRichData_local->getVertexProperty<double>("protein_velocity")
    //         .reinterpretTo(*mesh);
  }
}

#ifdef MEM3DG_WITH_NETCDF
std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>, EigenVectorX1d,
           EigenVectorX3dr, double>
System::readTrajFile(std::string trajFile, int startingFrame) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;
  double initialTime;
  EigenVectorX3dr initialVelocity;
  EigenVectorX1d initialProteinDensity;

  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  fd.getNcFrame(startingFrame);
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      fd.getCoords(startingFrame), fd.getTopology(startingFrame));
  const EigenVectorX3dr refVertexPositions = fd.getRefCoords(startingFrame);
  refVpg =
      std::make_unique<gcs::VertexPositionGeometry>(*mesh, refVertexPositions);

  // Map continuation variables
  initialTime = fd.getTime(startingFrame);
  initialVelocity = fd.getVelocity(startingFrame);
  initialProteinDensity = fd.getProteinDensity(startingFrame);
  // F.toMatrix(vel_protein) = fd.getProteinVelocity(startingFrame);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg),
                         initialProteinDensity, initialVelocity, initialTime);
}
#endif

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshFile(std::string inputMesh) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::readManifoldSurfaceMesh(inputMesh);

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMeshFile(std::string inputMesh, std::string referenceMesh) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::ManifoldSurfaceMesh> refMesh;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::readManifoldSurfaceMesh(inputMesh);
  std::tie(refMesh, refVpg) = gcs::readManifoldSurfaceMesh(referenceMesh);
  refVpg = refVpg->reinterpretTo(*mesh);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMatrices(EigenVectorX3sr &faceVertexMatrix,
                     EigenVectorX3dr &vertexPositionMatrix) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      vertexPositionMatrix, faceVertexMatrix);

  return std::make_tuple(std::move(mesh), std::move(vpg));
}

std::tuple<std::unique_ptr<gcs::ManifoldSurfaceMesh>,
           std::unique_ptr<gcs::VertexPositionGeometry>,
           std::unique_ptr<gcs::VertexPositionGeometry>>
System::readMatrices(EigenVectorX3sr &faceVertexMatrix,
                     EigenVectorX3dr &vertexPositionMatrix,
                     EigenVectorX3dr &refVertexPositionMatrix) {

  // Declare pointers to mesh / geometry objects
  std::unique_ptr<gcs::ManifoldSurfaceMesh> mesh;
  std::unique_ptr<gcs::VertexPositionGeometry> vpg;
  std::unique_ptr<gcs::VertexPositionGeometry> refVpg;

  // Load input mesh and geometry
  std::tie(mesh, vpg) = gcs::makeManifoldSurfaceMeshAndGeometry(
      vertexPositionMatrix, faceVertexMatrix);
  refVpg = std::make_unique<gcs::VertexPositionGeometry>(
      *mesh, refVertexPositionMatrix);

  return std::make_tuple(std::move(mesh), std::move(vpg), std::move(refVpg));
}

} // namespace solver
} // namespace mem3dg
