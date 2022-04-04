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

#pragma once

// #include <cassert>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include <Eigen/Core>

#include <math.h>
#include <vector>

#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include "mem3dg/constants.h"
#include "mem3dg/macros.h"
#include "mem3dg/mesh_io.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/type_utilities.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

namespace solver {

struct Forces {
  /// Cached mesh of interest
  gcs::ManifoldSurfaceMesh &mesh;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;

  /// Cached bending force
  gcs::VertexData<double> bendingForce;
  /// Cached deviatoric force
  gcs::VertexData<double> deviatoricForce;
  /// Cached osmotic force
  gcs::VertexData<double> osmoticForce;
  /// Cached tension-induced capillary force
  gcs::VertexData<double> capillaryForce;
  /// Cached interfacial line tension force
  gcs::VertexData<double> lineCapillaryForce;
  /// Cached adsorption induced surface force
  gcs::VertexData<double> adsorptionForce;
  /// Cached aggregation induced surface force
  gcs::VertexData<double> aggregationForce;
  /// Cached externally-applied force
  gcs::VertexData<double> externalForce;
  /// Cached self avoidance force
  gcs::VertexData<double> selfAvoidanceForce;
  /// Cached mechanical force
  gcs::VertexData<double> mechanicalForce;
  /// Cached osmotic pressure
  double osmoticPressure;
  /// Cached Surface tension
  double surfaceTension;

  /// Cached bending force
  gcs::VertexData<gc::Vector3> bendingForceVec;
  /// Cached bending force areaGrad component
  gcs::VertexData<gc::Vector3> bendingForceVec_areaGrad;
  /// Cached bending force gaussVec component
  gcs::VertexData<gc::Vector3> bendingForceVec_gaussVec;
  /// Cached bending force schlafliVec component
  gcs::VertexData<gc::Vector3> bendingForceVec_schlafliVec;

  /// Cached deviatoric force
  gcs::VertexData<gc::Vector3> deviatoricForceVec;
  /// Cached deviatoric force
  gcs::VertexData<gc::Vector3> deviatoricForceVec_mean;
  /// Cached deviatoric force
  gcs::VertexData<gc::Vector3> deviatoricForceVec_gauss;

  /// Cached tension-induced capillary force
  gcs::VertexData<gc::Vector3> capillaryForceVec;
  /// Cached osmotic force
  gcs::VertexData<gc::Vector3> osmoticForceVec;
  /// Cached Dirichlet energy driven force
  gcs::VertexData<gc::Vector3> lineCapillaryForceVec;
  /// Cached adsorption driven force
  gcs::VertexData<gc::Vector3> adsorptionForceVec;
  /// Cached aggregation driven force
  gcs::VertexData<gc::Vector3> aggregationForceVec;
  /// Cached external localized force
  gcs::VertexData<gc::Vector3> externalForceVec;
  /// Cached self avoidance force
  gcs::VertexData<gc::Vector3> selfAvoidanceForceVec;
  /// Cached mechanical force
  gcs::VertexData<gc::Vector3> mechanicalForceVec;

  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForceVec;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForceVec;

  /// Cached local stretching forces (in-plane regularization)
  gcs::VertexData<gc::Vector3> regularizationForce;

  /// Cached interior penalty chemical potential
  gcs::VertexData<double> interiorPenaltyPotential;
  /// Cached bending related chemical potential
  gcs::VertexData<double> bendingPotential;
  /// Cached deviatoric related chemical potential
  gcs::VertexData<double> deviatoricPotential;
  /// Cached adsorption related chemical potential
  gcs::VertexData<double> adsorptionPotential;
  /// Cached dirichlet energy related chemical potential
  gcs::VertexData<double> diffusionPotential;
  /// Cached aggregation related chemical potential
  gcs::VertexData<double> aggregationPotential;
  /// Cached chemical potential
  gcs::VertexData<double> chemicalPotential;

  /// force mask
  gcs::VertexData<gc::Vector3> forceMask;
  /// protein mask
  gcs::VertexData<double> proteinMask;

  Forces(gcs::ManifoldSurfaceMesh &mesh_, gcs::VertexPositionGeometry &vpg_)
      : mesh(mesh_), vpg(vpg_), mechanicalForce(mesh, 0),
        mechanicalForceVec(mesh, {0, 0, 0}), bendingForceVec(mesh, {0, 0, 0}),
        deviatoricForceVec(mesh, {0, 0, 0}),
        deviatoricForceVec_mean(mesh, {0, 0, 0}),
        deviatoricForceVec_gauss(mesh, {0, 0, 0}),
        bendingForceVec_areaGrad(mesh, {0, 0, 0}),
        bendingForceVec_gaussVec(mesh, {0, 0, 0}),
        bendingForceVec_schlafliVec(mesh, {0, 0, 0}),
        capillaryForceVec(mesh, {0, 0, 0}), osmoticForceVec(mesh, {0, 0, 0}),
        adsorptionForceVec(mesh, {0, 0, 0}),
        aggregationForceVec(mesh, {0, 0, 0}), externalForceVec(mesh, {0, 0, 0}),
        selfAvoidanceForceVec(mesh, {0, 0, 0}),
        lineCapillaryForceVec(mesh, {0, 0, 0}), bendingForce(mesh, 0),
        deviatoricForce(mesh, 0), capillaryForce(mesh, 0), surfaceTension(0),
        lineCapillaryForce(mesh, 0), adsorptionForce(mesh, 0),
        aggregationForce(mesh, 0), externalForce(mesh, 0),
        selfAvoidanceForce(mesh, 0), osmoticForce(mesh, 0), osmoticPressure(0),
        regularizationForce(mesh, {0, 0, 0}),
        stochasticForceVec(mesh, {0, 0, 0}), dampingForceVec(mesh, {0, 0, 0}),
        interiorPenaltyPotential(mesh, 0), bendingPotential(mesh, 0),
        deviatoricPotential(mesh, 0), adsorptionPotential(mesh, 0),
        aggregationPotential(mesh, 0), diffusionPotential(mesh, 0),
        chemicalPotential(mesh, 0), forceMask(mesh, {1.0, 1.0, 1.0}),
        proteinMask(mesh, 1) {}

  ~Forces() {}

  // ==========================================================
  // =============      Data interop helpers    ===============
  // ==========================================================

  /**
   * @brief Find the corresponding scalar vector (vertexData) by rowwise
   * (vertexwise) projecting matrix (vector vertexData) onto angle-weighted
   * normal vector
   */
  gcs::VertexData<double>
  ontoNormal(const gcs::VertexData<gc::Vector3> &vector) const {
    gcs::VertexData<double> vertexData(
        mesh, rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                gc::EigenMap<double, 3>(vpg.vertexNormals)));
    return vertexData;
  }
  gcs::VertexData<double>
  ontoNormal(const gcs::VertexData<gc::Vector3> &&vector) const {
    gcs::VertexData<double> vertexData(
        mesh, rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                gc::EigenMap<double, 3>(vpg.vertexNormals)));
    return vertexData;
  }

  EigenVectorX1d ontoNormal(const EigenVectorX3dr &vector) const {
    return rowwiseDotProduct(vector,
                             gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  EigenVectorX1d ontoNormal(const EigenVectorX3dr &&vector) const {
    return rowwiseDotProduct(vector,
                             gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  double ontoNormal(const gc::Vector3 &vector, const gc::Vertex &v) const {
    return gc::dot(vector, vpg.vertexNormals[v]);
  }

  double ontoNormal(const gc::Vector3 &vector, const std::size_t i) const {
    return gc::dot(vector, vpg.vertexNormals[i]);
  }

  double ontoNormal(const gc::Vector3 &&vector, gc::Vertex &v) const {
    return gc::dot(vector, vpg.vertexNormals[v]);
  }

  /**
   * @brief Find the corresponding matrix (vector vertexData) by rowwise
   * (vertexwise) appending angle-weighted normal vector to the scalar vector
   * (vertexData)
   */
  gcs::VertexData<gc::Vector3> addNormal(gcs::VertexData<double> &vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = rowwiseScalarProduct(
        vector.raw(), gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }
  gcs::VertexData<gc::Vector3> addNormal(gcs::VertexData<double> &&vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData) = rowwiseScalarProduct(
        vector.raw(), gc::EigenMap<double, 3>(vpg.vertexNormals));
    return vertexData;
  }

  EigenVectorX3dr addNormal(EigenVectorX1d &vector) {
    return rowwiseScalarProduct(vector,
                                gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  EigenVectorX3dr addNormal(EigenVectorX1d &&vector) {
    return rowwiseScalarProduct(vector,
                                gc::EigenMap<double, 3>(vpg.vertexNormals));
  }

  gc::Vector3 addNormal(double &vector, gc::Vertex &v) {
    return vector * vpg.vertexNormals[v];
  }
  gc::Vector3 addNormal(double &&vector, gc::Vertex &v) {
    return vector * vpg.vertexNormals[v];
  }

  /**
   * @brief Project the vector onto tangent plane by removing the
   * angle-weighted normal component
   */
  EigenVectorX3dr toTangent(gcs::VertexData<gc::Vector3> &vector) {
    return gc::EigenMap<double, 3>(vector) -
           rowwiseScalarProduct(
               rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                 gc::EigenMap<double, 3>(vpg.vertexNormals)),
               gc::EigenMap<double, 3>(vpg.vertexNormals));
  }
  EigenVectorX3dr toTangent(gcs::VertexData<gc::Vector3> &&vector) {
    return gc::EigenMap<double, 3>(vector) -
           rowwiseScalarProduct(
               rowwiseDotProduct(gc::EigenMap<double, 3>(vector),
                                 gc::EigenMap<double, 3>(vpg.vertexNormals)),
               gc::EigenMap<double, 3>(vpg.vertexNormals));
  }

  gc::Vector3 toTangent(gc::Vector3 &vector, gc::Vertex &v) {
    return vector -
           gc::dot(vector, vpg.vertexNormals[v]) * vpg.vertexNormals[v];
  }
  gc::Vector3 toTangent(gc::Vector3 &&vector, gc::Vertex &v) {
    return vector -
           gc::dot(vector, vpg.vertexNormals[v]) * vpg.vertexNormals[v];
  }

  /**
   * @brief Find the masked force
   */
  gcs::VertexData<gc::Vector3> maskForce(gcs::VertexData<gc::Vector3> &vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData).array() =
        gc::EigenMap<double, 3>(vector).array() *
        gc::EigenMap<double, 3>(forceMask).array();
    return vertexData;
  }

  gcs::VertexData<gc::Vector3>
  maskForce(gcs::VertexData<gc::Vector3> &&vector) {
    gcs::VertexData<gc::Vector3> vertexData(mesh);
    gc::EigenMap<double, 3>(vertexData).array() =
        gc::EigenMap<double, 3>(vector).array() *
        gc::EigenMap<double, 3>(forceMask).array();
    return vertexData;
  }

  EigenVectorX3dr maskForce(const EigenVectorX3dr &vector) const {
    return vector.array() * gc::EigenMap<double, 3>(forceMask).array();
  }

  EigenVectorX3dr maskForce(const EigenVectorX3dr &&vector) const {
    return vector.array() * gc::EigenMap<double, 3>(forceMask).array();
  }

  gc::Vector3 maskForce(const gc::Vector3 &vector, const gc::Vertex &v) const {
    return gc::Vector3{vector.x * forceMask[v].x, vector.y * forceMask[v].y,
                       vector.z * forceMask[v].z};
  }

  gc::Vector3 maskForce(const gc::Vector3 &vector, const std::size_t i) const {
    return gc::Vector3{vector.x * forceMask[i].x, vector.y * forceMask[i].y,
                       vector.z * forceMask[i].z};
  }

  gc::Vector3 maskForce(const gc::Vector3 &&vector, const gc::Vertex &v) const {
    return gc::Vector3{vector.x * forceMask[v].x, vector.y * forceMask[v].y,
                       vector.z * forceMask[v].z};
  }

  /**
   * @brief Find the masked chemical potential
   */
  gcs::VertexData<double> maskProtein(gcs::VertexData<double> &potential) {
    gcs::VertexData<double> vertexData(mesh);
    toMatrix(vertexData).array() =
        toMatrix(potential).array() * toMatrix(proteinMask).array();
    return vertexData;
  }

  gcs::VertexData<double> maskProtein(gcs::VertexData<double> &&potential) {
    gcs::VertexData<double> vertexData(mesh);
    toMatrix(vertexData).array() =
        toMatrix(potential).array() * toMatrix(proteinMask).array();
    return vertexData;
  }

  EigenVectorX1d maskProtein(EigenVectorX1d &potential) {
    return potential.array() * toMatrix(proteinMask).array();
  }

  EigenVectorX1d maskProtein(EigenVectorX1d &&potential) {
    return potential.array() * toMatrix(proteinMask).array();
  }

  double maskProtein(double &potential, gc::Vertex &v) {
    return potential * proteinMask[v];
  }

  double maskProtein(double &&potential, gc::Vertex &v) {
    return potential * proteinMask[v];
  }
};

} // namespace solver
} // namespace mem3dg
