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

#include <cstddef>
#include <geometrycentral/surface/geometry.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "Eigen/src/Core/util/Constants.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

#include "mem3dg/meshops.h"
#include "mem3dg/solver/mutable_trajfile.h"
#include "mem3dg/solver/trajfile.h"

#include <csignal>
#include <stdexcept>

namespace mem3dg {
namespace solver {
namespace integrator {

// ==========================================================
// =============        Integrator             ==============
// ==========================================================
class DLL_PUBLIC Integrator {
public:
  // variables (read-only)
  /// initial maximum force
  double initialMaximumForce;
  /// ratio of time step to the squared mesh size
  double dt_size2_ratio;
  /// Flag of success of the simulation
  bool SUCCESS = true;
  /// Flag for terminating the simulation
  bool EXIT = false;
  /// time step
  double timeStep;
  /// frame index
  size_t frame = 0;
  /// is continuing simulation
  bool isContinuation;
  /// whether disable the use of integrate()
  bool ifDisableIntegrate;
  /// System object to be integrated
  System &system;
  /// TrajFile
#ifdef MEM3DG_WITH_NETCDF
  MutableTrajFile mutableTrajFile;
  bool isMutableNetcdfFileCreated = false;
#endif

  // key parameters (read/write)
  /// characteristic time step
  double characteristicTimeStep;
  // total simulation time
  double totalTime = std::numeric_limits<double>::max();
  /// period of saving output data
  double savePeriod;
  /// tolerance for termination
  double tolerance;
  /// option to scale time step according to mesh size
  bool ifAdaptiveStep = false;
  /// path to the output directory
  std::string outputDirectory;
  /// option to use backtracking line search algorithm
  bool isBacktrack = true;
  /// backtracking coefficient
  double rho = 0.7;
  /// Wolfe condition parameter
  double c1 = 0.001;

  // defaulted parameters (read/write)
  /// if just save geometry .ply file
  bool ifJustGeometryPly = false;
  /// if output netcdf traj file
  bool ifOutputTrajFile = false;
  /// if output .ply file
  bool ifOutputMeshFile = false;
  /// if print to console
  bool ifPrintToConsole = false;
  /// name of the trajectory file
  std::string trajFileName = "traj.nc";

  // ==========================================================
  // =============        Constructor            ==============
  // ==========================================================
  /**
   * @brief Construct a new integrator object
   * @param f, System object to be integrated
   * @param dt_, characteristic time step
   * @param total_time, total simulation time
   * @param tSave, period of saving output data
   * @param tolerance, tolerance for termination
   * @param outputDir, path to the output directory
   */
  Integrator(System &system_, double characteristicTimeStep_, double totalTime_,
             double savePeriod_, double tolerance_,
             std::string outputDirectory_, std::size_t frame_ = 0)
      : Integrator(system_, characteristicTimeStep_, tolerance_,
                   outputDirectory_) {
    totalTime = totalTime_;
    frame = frame_;
    savePeriod = savePeriod_;

    ifDisableIntegrate = false;
    ifPrintToConsole = false;
    isContinuation = (frame != 0);
  }

  Integrator(System &system_, double characteristicTimeStep_, double tolerance_,
             std::string outputDirectory_)
      : system(system_), characteristicTimeStep(characteristicTimeStep_),
        tolerance(tolerance_), timeStep(characteristicTimeStep_),
        outputDirectory(outputDirectory_) {
    ifDisableIntegrate = true;
    ifPrintToConsole = true;
    // Initialize the timestep-meshsize ratio
    dt_size2_ratio =
        characteristicTimeStep /
        std::pow(system.geometry.vpg->edgeLengths.raw().minCoeff(), 2);

    // Initialize the initial maxForce
    system.computeConservativeForcing();
    system.addNonconservativeForcing(timeStep);
    initialMaximumForce =
        system.parameters.variation.isShapeVariation
            ? system.forces.mechanicalForce.raw().cwiseAbs().maxCoeff()
            : system.forces.chemicalPotential.raw().cwiseAbs().maxCoeff();
  }

  /**
   * @brief Destroy the Integrator
   *
   */
  virtual ~Integrator() = default;

  // ==========================================================
  // =================   Template functions    ================
  // ==========================================================
  virtual bool integrate() = 0;
  virtual void march() = 0;
  virtual void status() = 0;
  virtual void checkParameters() = 0;

  // ==========================================================
  // =================     Output Data         ================
  // ==========================================================
  /**
   * @brief Save trajectory, mesh and print to console
   */
  void saveData(bool ifOutputTrajFile, bool ifOutputMeshFile,
                bool ifPrintToConsole);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Initialize netcdf trajectory file
   */
  void createMutableNetcdfFile(bool isContinue);

  /**
   * @brief close netcdf trajectory file
   */
  void closeMutableNetcdfFile();

  /**
   * @brief Save data to netcdf trajectory file
   */
  void saveMutableNetcdfData();

#endif

  // ==========================================================
  // =============     Helper functions          ==============
  // ==========================================================

  /**
   * @brief Backtracking algorithm that dynamically adjust step size based on
   * energy evaluation
   * @param positionDirection, direction of shape, most likely some function of
   * gradient
   * @param chemicalDirection, direction of protein density, most likely some
   * function of gradient
   * @param rho, discount factor
   * @param c1, constant for Wolfe condition, between 0 to 1, usually ~ 1e-4
   * @return alpha, line search step size
   */
  double backtrack(Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection,
                   double rho = 0.7, double c1 = 0.001);

  /**
   * @brief Backtracking algorithm that dynamically adjust step size based on
   * energy evaluation
   * @param positionDirection, direction of shape, most likely some function of
   * gradient
   * @param rho, discount factor
   * @param c1, constant for Wolfe condition, between 0 to 1, usually ~ 1e-4
   * @return alpha, line search step size
   */
  double mechanicalBacktrack(
      Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
      double rho = 0.7, double c1 = 0.001);

  /**
   * @brief Backtracking algorithm that dynamically adjust step size based on
   * energy evaluation
   * @param chemicalDirection, direction of protein density, most likely some
   * function of gradient
   * @param rho, discount factor
   * @param c1, constant for Wolfe condition, between 0 to 1, usually ~ 1e-4
   * @return alpha, line search step size
   */
  double
  chemicalBacktrack(Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection,
                    double rho = 0.7, double c1 = 0);

  /**
   * @brief get adaptive characteristic time step
   * @return characteristic time step
   */
  double getAdaptiveCharacteristicTimeStep();

  // /**
  //  * @brief fix point iteration implicit solve
  //  * @param FlowMap functional of position and protein flow map
  //  * @param h time step
  //  * @param tolearance convergence tolerance relative to timestep
  //  */
  // std::tuple<EigenVectorX3dr, EigenVectorX1d>
  // fixPointIteration(std::function<std::tuple<EigenVectorX3dr,
  // EigenVectorX1d>(
  //                       EigenVectorX3dr &, EigenVectorX1d &, double)>
  //                       flowMap,
  //                   const double h, const double tolereance);
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
