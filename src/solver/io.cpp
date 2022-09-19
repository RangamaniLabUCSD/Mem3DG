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
    gcs::writeSurfaceMesh(*geometry.mesh, *geometry.vpg, PathToSave);
  } else {
    gcs::RichSurfaceMeshData richData(*geometry.mesh);
    richData.addMeshConnectivity();
    richData.addGeometry(*geometry.vpg);

    // write protein distribution
    richData.addVertexProperty("proteinDensity", proteinDensity);
    richData.addVertexProperty("velocity", forces.ontoNormal(velocity));

    // write bool
    gcs::VertexData<double> msk(*geometry.mesh);
    msk.fromVector(toMatrix(forces.forceMask).rowwise().sum());
    richData.addVertexProperty("forceMask", msk);
    richData.addVertexProperty("proteinMask", forces.proteinMask);
    gcs::VertexData<double> tkr(*geometry.mesh);
    tkr.fromVector(geometry.notableVertex.raw().cast<double>());
    richData.addVertexProperty("notableVertex", tkr);
    // gcs::VertexData<int> mutMkr(*geometry.mesh);
    // mutMkr.fromVector(mutationMarker.raw().cast<int>());
    // richData.addVertexProperty("smoothing_mask", mutMkr);

    // write geometry
    gcs::VertexData<double> meanCurv(*geometry.mesh);
    meanCurv.fromVector(geometry.vpg->vertexMeanCurvatures.raw().array() /
                        geometry.vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("meanCurvature", meanCurv);
    gcs::VertexData<double> gaussCurv(*geometry.mesh);
    gaussCurv.fromVector(geometry.vpg->vertexGaussianCurvatures.raw().array() /
                         geometry.vpg->vertexDualAreas.raw().array());
    richData.addVertexProperty("gaussianCurvature", gaussCurv);
    richData.addVertexProperty("spontaneousCurvature", H0);
    richData.addVertexProperty("dualArea", geometry.vpg->vertexDualAreas);

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
  if (geometry.mesh->nFaces() != ptrMesh_local->nFaces() ||
      geometry.mesh->nVertices() != ptrMesh_local->nVertices()) {
    throw std::logic_error(
        "Topology for continuation parameters mapping is not consistent!");
  } else {
    // Map continuation variables
    time = std::stod(sliceString(plyFile, "t", "_"));
    energy.time = time;
    // frame = std::stod(sliceString(plyFile, "f", "_"));
    proteinDensity =
        ptrRichData_local->getVertexProperty<double>("protein_density")
            .reinterpretTo(*geometry.mesh);
    // vel_protein =
    //     ptrRichData_local->getVertexProperty<double>("protein_velocity")
    //         .reinterpretTo(*geometry.mesh);
  }
}

#ifdef MEM3DG_WITH_NETCDF
// std::tuple<Geometry &&, EigenVectorX1d &, EigenVectorX3dr &, double>
std::tuple<EigenVectorX1d, EigenVectorX3dr, double>
System::readTrajFile(std::string trajFile, int startingFrame) {
  // Geometry geometry_here(trajFile, startingFrame);
  // std::unique_ptr<Geometry> geometry_here =
  //     std::make_unique<Geometry>(trajFile, startingFrame);
  MutableTrajFile fd = MutableTrajFile::openReadOnly(trajFile);
  // Map continuation variables
  double time = fd.getTime(startingFrame);
  EigenVectorX3dr initialVelocity = fd.getVelocity(startingFrame);
  EigenVectorX1d initialProteinDensity = fd.getProteinDensity(startingFrame);
  // F.toMatrix(vel_protein) = fd.getProteinVelocity(startingFrame);

  return std::make_tuple(initialProteinDensity, initialVelocity, time);
  // return std::forward_as_tuple(Geometry(trajFile, startingFrame),
  // fd.getProteinDensity(startingFrame),  fd.getVelocity(startingFrame), time);
}
#endif

} // namespace solver
} // namespace mem3dg
