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
protected:
  /// Energy of the last time step
  Energy previousE;
  /// last time saving the data
  double lastSave;
  /// last time updating geodesics
  double lastUpdateGeodesics;
  /// last time processing mesh
  double lastProcessMesh;
  /// regularization force to the system
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce;
  /// numerical dissipative particle dynamics force to the system
  Eigen::Matrix<double, Eigen::Dynamic, 1> DPDForce;
  /// Starting time of the simulation
  double init_time;
  /// Flag of success of the simulation
  bool SUCCESS;
  /// Flag for terminating the simulation
  bool EXIT;
  /// Frame index of the trajectory output
  std::size_t frame;
  /// Normalized area difference to reference mesh
  double dArea;
  /// Normalized volume/osmotic pressure difference
  double dVP;
  /// ratio of time step to the squared mesh size
  double dt_size2_ratio;
  /// initial maximum force
  double maxForce;
  /// TrajFile
#ifdef MEM3DG_WITH_NETCDF
  TrajFile trajFile;
  MutableTrajFile mutableTrajFile;
#endif

public:
  /// System object to be integrated
  System &f;
  /// time step
  double dt;
  // total simulation time
  double total_time;
  /// period of saving output data
  double tSave;
  /// tolerance for termination
  double tol;
  /// path to the output directory
  std::string outputDir;

  /// period of saving output data
  double tUpdateGeodesics;
  /// period of saving output data
  double tProcessMesh;
  /// name of the trajectory file
  std::string trajFileName = "traj.nc";
  /// option to scale time step according to mesh size
  bool isAdaptiveStep = true;
  /// verbosity level of integrator
  size_t verbosity = 3;
  /// just save geometry .ply file
  bool isJustGeometryPly = false;

  // ==========================================================
  // =============        Constructor            ==============
  // ==========================================================
  /**
   * @brief Construct a new integrator object
   * @param f, System object to be integrated
   * @param dt, time step
   * @param total_time, total simulation time
   * @param tSave, period of saving output data
   * @param tolerance, tolerance for termination
   * @param outputDir, path to the output directory
   */
  Integrator(System &f_, double dt_, double total_time_, double tSave_,
             double tolerance_, std::string outputDir_)
      : f(f_), previousE(f_.energy), dt(dt_), total_time(total_time_),
        tSave(tSave_), tol(tolerance_), tUpdateGeodesics(total_time),
        tProcessMesh(total_time), outputDir(outputDir_), init_time(f_.time),
        SUCCESS(true), EXIT(false), frame(0), lastUpdateGeodesics(f_.time),
        lastProcessMesh(f_.time), lastSave(f_.time) {

    // Initialize the timestep-meshsize ratio
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();

    // Initialize the initial maxForce
    getForces();
    maxForce = f.parameters.variation.isShapeVariation
                   ? toMatrix(f.forces.mechanicalForce).cwiseAbs().maxCoeff()
                   : f.forces.chemicalPotential.raw().cwiseAbs().maxCoeff();

    // Initialize geometry constraints
    dArea = 1e10;
    dVP = 1e10;

    // Initialize system summarized forces
    regularizationForce.resize(f.mesh->nVertices(), 3);
    DPDForce.resize(f.mesh->nVertices(), 1);
    regularizationForce.setZero();
    DPDForce.setZero();
  }

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
  void saveData();

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Initialize netcdf traj file
   */
  void createNetcdfFile();
  /**
   * @brief Initialize netcdf traj file
   */
  void createMutableNetcdfFile();
  /**
   * @brief Save data to netcdf traj file
   */
  void saveNetcdfData();
  /**
   * @brief Save data to netcdf traj file
   */
  void saveMutableNetcdfData();

#endif

  /**
   * @brief Mark the file name
   *
   * @param dirPath path of the directory
   * @param file name of the file, for example in the form of "/traj.nc"
   * @param marker_str marker used to mark the file, such as marker = "_failed"
   * results in new file name of "/traj_failed.nc"
   */
  void markFileName(std::string marker_str);

  /**
   * @brief Save parameters to txt file
   * @return
   */
  void getParameterLog(std::string inputMesh);

  /**
   * @brief Save status log to txt file
   * @return
   */
  void getStatusLog(std::string nameOfFile, std::size_t frame, double areaError,
                    double volumeError, double bendingError, double faceError,
                    std::string inputMesh);

  // ==========================================================
  // =============     Helper functions          ==============
  // ==========================================================
  /**
   * @brief Thresholding when adopting reduced volume parametrization
   * @param EXIT, reference to the exit flag
   * @param isAugmentedLagrangian, whether using augmented lagrangian method
   * @param dArea, normalized area difference
   * @param dVolume, normalized volume difference
   * @param ctol, exit criterion for constraint
   * @param increment, increment coefficient of penalty when using incremental
   * penalty method
   * @return
   */
  void reducedVolumeThreshold(bool &EXIT, const bool isAugmentedLagrangian,
                              const double dArea, const double dVolume,
                              const double ctol, double increment);
  /**
   * @brief Thresholding when adopting ambient pressure constraint
   * @param EXIT, reference to the exit flag
   * @param isAugmentedLagrangian, whether using augmented lagrangian method
   * @param dArea, normalized area difference
   * @param ctol, exit criterion for constraint
   * @param increment, increment coefficient of penalty when using incremental
   * penalty method
   * @return
   */
  void pressureConstraintThreshold(bool &EXIT, const bool isAugmentedLagrangian,
                                   const double dArea, const double ctol,
                                   double increment);

  /**
   * @brief Backtracking algorithm that dynamically adjust step size based on
   * energy evaluation
   * @param potentialEnergy_pre, previous energy evaluation
   * @param positionDirection, direction of shape, most likely some function of
   * gradient
   * @param chemicalDirection, direction of protein density, most likely some
   * function of gradient
   * @param rho, discount factor
   * @param c1, constant for Wolfe condtion, between 0 to 1, usually ~ 1e-4
   * @return alpha, line search step size
   */
  double backtrack(const double energy_pre,
                   Eigen::Matrix<double, Eigen::Dynamic, 3> &&positionDirection,
                   Eigen::Matrix<double, Eigen::Dynamic, 1> &chemicalDirection,
                   double rho = 0.99, double c1 = 0.0001);

  /**
   * @brief Summerize forces into 3 categories: physcialPressure, DPDPressure
   * and regularizationForce. Note that the forces has been removed rigid body
   * mode and masked for integration
   * @return
   */
  void getForces();

  /**
   * @brief Check finiteness of simulation states and backtrack for error in
   * specific component
   * @return
   */
  void finitenessErrorBacktrack();

  /**
   * @brief Backtrack the line search failure by inspecting specific
   * energy-force relation
   * @return
   */
  void lineSearchErrorBacktrack(const double &alpha,
                                const EigenVectorX3dr &initial_pos,
                                const EigenVectorX1d &init_proteinDensity,
                                bool runAll = false);
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
