#pragma once

#include <cassert>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/heat_method_distance.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

#include "ddgsolver/util.h"
#include "ddgsolver/meshops.h"
#include "ddgsolver/macros.h"

namespace ddgsolver {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

struct DLL_PUBLIC Parameters {
  /// Bending modulus
  double Kb;   
  /// Spontaneous curvature
  double H0;    
  /// Local stretching modulus
  double Ksl;   
  /// Global stretching modulus
  double Ksg;   
  /// Edge spring constant
  double Kse; 
  /// Volume regularization
  double Kv;   
  /// Dissipation coefficient
  double gamma; 
  /// Reduced volume
  double Vt;    
  /// Boltzmann constant*Temperature
  double kt;  
  /// Noise
  double sigma; 
  /// index of node with applied external force 
  int ptInd;  
  /// Magnitude of external force 
  double extF;   
  /// level of concentration of the external force
  double conc;   
};

class DLL_PUBLIC Force : public Parameters {
public:
  /// Cached mesh of interest
  gcs::HalfedgeMesh &mesh;
  /// Embedding and other geometric details
  gcs::VertexPositionGeometry &vpg;
  /// Cached bending forces
  gcs::VertexData<gc::Vector3> bendingForces;
  /// Cached stretching forces
  gcs::VertexData<gc::Vector3> stretchingForces;
  /// Cached pressure induced forces
  gcs::VertexData<gc::Vector3> pressureForces;
  /// Cached damping forces
  gcs::VertexData<gc::Vector3> dampingForces;
  /// Cached stochastic forces
  gcs::VertexData<gc::Vector3> stochasticForces;
  /// Cached external forces
  gcs::VertexData<gc::Vector3> externalForces;

  /// Target area per face
  gcs::FaceData<double> initialFaceAreas;
  /// Target total face area
  double initialSurfaceArea = 0.0;
  /// surface area
  double surfaceArea = 0.0;
  /// Target length per edge
  gcs::EdgeData<double> targetEdgeLength;
  /// Maximal volume
  double maxVolume = 0.0;
  /// Volume
  double volume = 0.0;
  /// Cached vertex positions from the previous step
  gcs::VertexData<gc::Vector3> pastPositions;
  /// Cached vertex velocity by finite differencing past and current position
  gcs::VertexData<gc::Vector3> vel;
  // Mean curvature of the mesh
  gcs::VertexData<double> H;
  // Spontaneous curvature of the mesh
  gcs::VertexData<double> H0;
  /// Random number engine
  pcg32 rng;
  std::normal_distribution<double> normal_dist;
  /// magnitude of externally applied force
  Eigen::Matrix<double, Eigen::Dynamic, 1> appliedForceMagnitude;

  /**
   * @brief Construct a new Force object
   *
   * @param mesh_         Mesh connectivity
   * @param vpg_          Embedding and geometry information
   * @param time_step_    Numerical timestep
   */

  Force(gcs::HalfedgeMesh &mesh_, gcs::VertexPositionGeometry &vpg_, Parameters &p)
      : mesh(mesh_), vpg(vpg_), bendingForces(mesh_, {0, 0, 0}), Parameters(p),
        stretchingForces(mesh_, {0, 0, 0}), dampingForces(mesh_, {0, 0, 0}),
        pressureForces(mesh_, {0, 0, 0}), stochasticForces(mesh_, {0, 0, 0}),
        externalForces(mesh_, {0, 0, 0}), vel(mesh_, { 0, 0, 0 }) {

    // Initialize RNG
    pcg_extras::seed_seq_from<std::random_device> seed_source;
    rng = pcg32(seed_source);

    // GC computed properties
    vpg.requireFaceNormals();
    vpg.requireVertexDualAreas();
    vpg.requireEdgeDihedralAngles();
    vpg.requireCornerAngles();
    vpg.requireFaceAreas();
    vpg.requireVertexIndices();
    vpg.requireVertexGaussianCurvatures();
    vpg.requireFaceIndices();
    vpg.requireEdgeLengths();
    vpg.requireVertexNormals();


    // Initialize face areas
    initialFaceAreas = vpg.faceAreas;
    auto faceAreas_e = EigenMap(initialFaceAreas);
    initialSurfaceArea = faceAreas_e.sum();

    // Initialize edge length
    targetEdgeLength = vpg.edgeLengths;

    // Initialize maximal volume
    double pi = 2 * std::acos(0.0);
    maxVolume = std::pow(initialSurfaceArea, 1.5) / std::pow(pi, 0.5) / 6;
    //for (gcs::Face f : mesh.faces()) {
    //  maxVolume += signedVolumeFromFace(f, vpg);
    //}

    // Initialize the vertex position of the last iteration
    pastPositions = vpg.inputVertexPositions;

    // Initialize mean curvature
    H.fill(0.0);
    // Initialize the magnitude of externally applied force
    gcs::VertexData<double> geodesicDistanceFromAppliedForce 
      = heatMethodDistance(vpg, mesh.vertex(ptInd));
    auto dist_e = geodesicDistanceFromAppliedForce.toMappedVector();
    double stdDev = dist_e.maxCoeff()/conc;
    appliedForceMagnitude = extF / (stdDev * pow(pi * 2, 0.5))
      * (-dist_e.array() * dist_e.array()
        / (2 * stdDev * stdDev)).exp();
  }


  /**
   * @brief Destroy the Force object
   *
   * Explicitly unrequire values required by the constructor. In case, there
   * is another pointer to the HalfEdgeMesh and VertexPositionGeometry
   * elsewhere, calculation of dependent quantities should be respected.
   */
  ~Force() {
    vpg.unrequireFaceNormals();
    vpg.unrequireEdgeDihedralAngles();
    vpg.unrequireVertexDualAreas();
    vpg.unrequireFaceAreas();
    vpg.unrequireVertexIndices();
    vpg.unrequireVertexGaussianCurvatures();
    vpg.unrequireFaceIndices();
    vpg.unrequireEdgeLengths();
    vpg.unrequireVertexNormals();
  }

  void getBendingForces();

  void getStretchingForces();

  void getPressureForces();

  void getConservativeForces();

  void getDPDForces();

  void getExternalForces();


  /**
   * @brief Get velocity from the position of the last iteration
   *
   * @param timeStep
   */
  void getVelocityFromPastPosition(double dt);

  /**
   * @brief Update the vertex position and recompute cached values
   *
   */
  void update_Vertex_positions() { vpg.refreshQuantities(); }

  void pcg_test();
};
} // end namespace ddgsolver
