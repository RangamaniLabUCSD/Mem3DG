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
#include <Eigen/SparseLU>

#include <pcg_random.hpp>
#include <random>

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
#include "mem3dg/solver/forces.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/type_utilities.h"

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace mem3dg {

namespace solver {

struct Parameters {
  struct Bending {
    /// membrane thickness
    double D = 0;
    /// unity modulus
    double alpha = 0;
    /// preferred area difference
    double dA0 = 0;

    /// Deviatoric modulus
    double Kd = 0;
    /// Constant of deviatoric modulus vs protein density
    double Kdc = 0;

    /// Bending modulus
    double Kb = 0;
    /// Constant of bending modulus vs protein density
    double Kbc = 0;
    /// Constant of Spontaneous curvature vs protein density
    double H0c = 0;
    /// type of relation between H0 and protein density, "linear" or "hill"
    std::string relation = "linear";

    /**
     * @brief check parameter conflicts
     */
    DLL_PUBLIC void checkParameters();
  };

  struct Tension {
    /// tension(area) function form
    std::function<std::tuple<double, double>(double)> form = NULL;
  };

  struct Osmotic {
    /// pressure(volume) function form
    std::function<std::tuple<double, double>(double)> form = NULL;
  };

  struct Adsorption {
    /// binding energy per protein
    double epsilon = 0;
  };

  struct Aggregation {
    /// aggregation energy constant
    double chi = 0;
  };

  struct Entropy {
    /// entropy energy constant
    double xi = 0;
  };

  struct Dirichlet {
    /// Smooothing coefficients
    double eta = 0;
  };

  struct SelfAvoidance {
    /// limit distance
    double d = 0;
    /// penalty coefficient
    double mu = 0;
    // neighborhood layers
    std::size_t n = 1;
    // period factor of computation
    double p = 0;
  };

  struct External {
    /// form of external force
    std::function<EigenVectorX3dr(EigenVectorX3dr, EigenVectorX1d, double,
                                  EigenVectorX1d)>
        form = NULL;
  };

  struct DPD {
    /// Dissipation coefficient
    double gamma = 0;
  };

  struct Boundary {
    /// shape boundary condition: roller, pin, fixed, none
    std::string shapeBoundaryCondition = "none";
    /// protein boundary condition: pin
    std::string proteinBoundaryCondition = "none";

    /**
     * @brief check parameter conflicts
     */
    DLL_PUBLIC void checkParameters();
  };

  struct Variation {
    /// Whether or not consider protein binding
    bool isProteinVariation = false;
    /// Whether conserve protein mass;
    bool isProteinConservation = false;
    /// Whether or not consider shape evolution
    bool isShapeVariation = true;
    /// domain of shape variation
    double geodesicMask = -1;
    /// period of updating mask
    std::size_t updateMaskPeriod = std::numeric_limits<std::size_t>::max();

    /**
     * @brief check parameter conflicts
     */
    DLL_PUBLIC void checkParameters();
  };

  struct Point {
    /// prescription of center finding
    std::function<Eigen::Matrix<bool, Eigen::Dynamic, 1>(
        EigenVectorX3sr, EigenVectorX3dr, EigenVectorX1d)>
        prescribeNotableVertex = NULL;
    /// period of updating geodesic distance from notableVertex calculation
    std::size_t updateGeodesicsPeriod = std::numeric_limits<std::size_t>::max();
    /// period of updating notable vertex based functional
    /// prescribeNotableVertex
    std::size_t updateNotableVertexPeriod =
        std::numeric_limits<std::size_t>::max();
  };

  struct Protein {
    /// interior point parameter for protein density
    double proteinInteriorPenalty = 0; // 1e-6
    /// precription of protein density
    std::function<EigenVectorX1d(double, EigenVectorX1d, EigenVectorX1d)>
        prescribeProteinDensityDistribution = NULL;
    /// period of updating protein density distribution
    std::size_t updateProteinDensityDistributionPeriod =
        std::numeric_limits<std::size_t>::max();
  };

  struct Spring {
    /// triangle ratio constant
    double Kst = 0;
    /// Local stretching modulus
    double Ksl = 0;
    /// Edge spring constant
    double Kse = 0;
  };

  /// bending parameters
  Bending bending;
  /// surface tension parameters
  Tension tension;
  /// osmotic pressure parameters
  Osmotic osmotic;
  /// protein adsorption parameters
  Adsorption adsorption;
  /// protein aggregation parameters
  Aggregation aggregation;
  /// protein entropy parameters
  Entropy entropy;
  /// protein dirichlet energy parameters
  Dirichlet dirichlet;
  /// self avoidance energy parameters
  SelfAvoidance selfAvoidance;
  /// external force parameters
  External external;
  /// DPD parameters
  DPD dpd;
  /// boundary conditions
  Boundary boundary;
  /// variation
  Variation variation;
  /// reference point
  Point point;
  /// protein distribution
  Protein protein = Protein();
  /// mesh regularizer
  Spring spring;

  /// mobility constant
  double proteinMobility = 0;
  /// Temperature
  double temperature = 293;
  /// damping coefficient
  double damping = 0;

  /**
   * @brief check parameter conflicts
   */
  DLL_PUBLIC void checkParameters(bool hasBoundary, size_t nVertex);
};

} // namespace solver
} // namespace mem3dg
