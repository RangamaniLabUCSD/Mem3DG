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

#include <geometrycentral/surface/geometry.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>

#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

#include "mem3dg/solver/meshops.h"
#include "system.h"
#include <csignal>
#ifdef MEM3DG_WITH_NETCDF
#include "mem3dg/solver/trajfile.h"
#endif

namespace mem3dg {

// ==========================================================
// =============        Integrator             ==============
// ==========================================================
class DLL_PUBLIC Integrator {
public:
  /// System object to be integrated
  System &f;
  /// time step
  double dt;
  // total simulation time
  double total_time;
  /// period of saving output data
  double tSave;
  /// verbosity level of integrator
  size_t verbosity;
  /// tolerance for termination
  double tol;
  /// path to the output directory
  std::string outputDir;
  /// option to scale time step according to mesh size
  bool isAdaptiveStep;
  /// name of the trajectory file
  std::string trajFileName;
  /// regularization force to the system
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce;
  /// physical vertex pressure to the system
  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure;
  /// numerical dissipative particle dynamics force to the system
  Eigen::Matrix<double, Eigen::Dynamic, 3> DPDPressure;

  /// Starting time of the simulation
  double init_time;
  /// Flag of success of the simulation
  bool SUCCESS;
  /// Flag for terminating the simulation
  bool EXIT;
  /// Frame index of the trajectory output
  size_t frame;
  /// Normalized area difference to reference mesh
  double dArea;
  /// Normalized volume/osmotic pressure difference
  double dVP;
  /// ratio of time step to the squared mesh size
  double dt_size2_ratio;

#ifdef MEM3DG_WITH_NETCDF
  TrajFile fd;
#endif

  // ==========================================================
  // =============        Constructor            ==============
  // ==========================================================
  /**
   * @brief Construct a new integrator object
   * @param f, System object to be integrated
   * @param dt, time step
   * @param isAdaptiveStep, option to scale time step according to mesh size
   * @param total_time, total simulation time
   * @param tSave, period of saving output data
   * @param tolerance, tolerance for termination
   * @param outputDir, path to the output directory
   * @param trajFileName, name of the trajectory file
   * @param verbosity,   0: no output file 1: Traj.nc 2: All above + console
   * printing 3: All above + txt + .ply
   */
  Integrator(System &f_, double dt_, bool isAdaptiveStep_, double total_time_,
             double tSave_, double tolerance_, std::string outputDir_,
             std::string trajFileName_, size_t verbosity_)
      : f(f_), dt(dt_), isAdaptiveStep(isAdaptiveStep_),
        total_time(total_time_), tSave(tSave_), tol(tolerance_),
        verbosity(verbosity_), outputDir(outputDir_),
        trajFileName(trajFileName_), init_time(f.time), SUCCESS(true),
        EXIT(false), frame(0) {

    // Initialize the timestep-meshsize ratio
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();

    // Initialize system summarized forces
    regularizationForce.resize(f.mesh->nVertices(), 3);
    physicalPressure.resize(f.mesh->nVertices(), 3);
    DPDPressure.resize(f.mesh->nVertices(), 3);

    // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
    createNetcdfFile();
#endif

    // print to console
    std::cout << "Initialized integrator and the output directory is "
              << outputDir << std::endl;
  }

  // ==========================================================
  // =================     Output Data         ================
  // ==========================================================
  /**
   * @brief Initialize netcdf traj file
   */
  void createNetcdfFile();

  /**
   * @brief Save trajectory, mesh and print to console
   */
  void saveData();

  /**
   * @brief Save data to richData
   */
  void saveRichData();

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Save data to netcdf traj file
   * @param frame, frame index of netcdf traj file
   * @param fd, netcdf trajFile object
   */
  void saveNetcdfData(size_t &frame, TrajFile &fd);

#endif
  /**
   * @brief Mark the file name
   *
   * @param dirPath path of the directory
   * @param file name of the file, for example in the form of "/traj.nc"
   * @param marker marker used to mark the file, such as marker = "_failed"
   * results in new file name of "/traj_failed.nc"
   */
  void markFileName(std::string marker_str) {
    std::string dirPath = outputDir;

    const char *marker = marker_str.c_str();

    char *file = new char[trajFileName.size() + 1];
    std::copy(trajFileName.begin(), trajFileName.end(), file);
    file[trajFileName.size()] = '\0';

    char fileMarked[50], oldNC[150], newNC[150];

    // sprintf(fileMarked, "/traj_H_%d_VP_%d_failed.nc", int(H * 100),
    //         int(VP * 100));

    // split the extension and file name
    char *ext;
    ext = strchr(file, '.');

    // name fileMarked to be the file name
    strncpy(fileMarked, file, ext - file);

    // name fileMarked to be file name + the marker + extension
    strcat(fileMarked, marker);
    strcat(fileMarked, ext);

    // append the directory path and copy to oldNC and newNC
    strcpy(oldNC, dirPath.c_str());
    strcpy(newNC, dirPath.c_str());
    strcat(oldNC, file);
    strcat(newNC, fileMarked);

    // rename file
    rename(oldNC, newNC);
    delete[] file;
  }
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
   * @param rho, discount factor
   * @param c1, constant for Wolfe condtion, between 0 to 1, usually ~ 1e-4
   * @param EXIT, exit flag for integration loop
   * @param SUCCESS, test flag on simulation sucess
   * @param potentialEnergy_pre, previous energy evaluation
   * @param force, gradient of the energy
   * @param direction, direction, most likely some function of gradient
   * @return
   */
  void backtrack(double rho, double c1, bool &EXIT, bool &SUCCESS,
                 const double potentialEnergy_pre,
                 const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
                 const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction);

  /**
   * @brief Summerize forces into 3 categories: physcialPressure, DPDPressure
   * and regularizationForce. Note that the forces has been removed rigid body
   * mode and masked for integration
   * @return
   */
  void getForces();

  /**
   * @brief Check parameters for time integration
   *
   * @param integrator, name of time integration
   * @return
   */
  void checkParameters(std::string integrator);
};

// ==========================================================
// =============        Stormer Verlet         ==============
// ==========================================================
/**
 * @brief Stomer Verlet time Integration
 */
class DLL_PUBLIC StormerVerlet : public Integrator {
public:
  StormerVerlet(System &f_, double dt_, bool isAdaptiveStep_,
                double total_time_, double tSave_, double tolerance_,
                std::string outputDir_, std::string trajFileName_,
                size_t verbosity_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_) {}

  void integrate();
};

// ==========================================================
// =============        Velocity Verlet        +=============
// ==========================================================
/**
 * @brief Velocity Verlet time Integration
 */
class DLL_PUBLIC VelocityVerlet : public Integrator {
public:
  // initialize variables used in time integration
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure;
  Eigen::Matrix<double, Eigen::Dynamic, 3> newTotalPressure;

  VelocityVerlet(System &f_, double dt_, bool isAdaptiveStep_,
                 double total_time_, double tSave_, double tolerance_,
                 std::string outputDir_, std::string trajFileName_,
                 size_t verbosity_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_) {
    if (verbosity > 1) {
      std::cout << "Running Velocity Verlet integrator ..." << std::endl;
    }
    // check the validity of parameter
    checkParameters("velocity verlet");
  }
  void integrate();

  /**
   * @brief velocity Verlet stepper
   */
  void step();

  /**
   * @brief Velocity Verlet status computation and thresholding
   */
  void status();
};

// ==========================================================
// =============        Forward Euler          +=============
// ==========================================================
/**
 * @brief Euler (gradient descent) time Integration
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @return Success, if simulation is sucessful
 */
class DLL_PUBLIC Euler : public Integrator {
public:
  const bool isBacktrack;
  const double rho;
  const double c1;
  Euler(System &f_, double dt_, bool isAdaptiveStep_, double total_time_,
        double tSave_, double tolerance_, std::string outputDir_,
        std::string trajFileName_, size_t verbosity_, bool isBacktrack_,
        double rho_, double c1_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_) {
    if (verbosity > 1) {
      std::cout << "Running Forward Euler (steepest descent) propagator ..."
                << std::endl;
    }
    // check the validity of parameter
    checkParameters("euler");
  }
  bool integrate();

  /**
   * @brief Forward Euler stepper
   */
  void step();

  /**
   * @brief Forward Euler status computation and thresholding
   */
  void status();
};

// ==========================================================
// =============      Conjugate Gradient        =============
// ==========================================================
/**
 * @brief Conjugate Gradient propagator
 * @param ctol, tolerance for termination (contraints)
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @return Success, if simulation is sucessful
 */
class DLL_PUBLIC ConjugateGradient : public Integrator {
public:
  const bool isBacktrack;
  const double rho;
  const double c1;
  double ctol;
  const bool isAugmentedLagrangian;

  Eigen::Matrix<double, Eigen::Dynamic, 3> direction;
  double currentNormSq;
  double pastNormSq;

  ConjugateGradient(System &f_, double dt_, bool isAdaptiveStep_,
                    double total_time_, double tSave_, double tolerance_,
                    std::string outputDir_, std::string trajFileName_,
                    size_t verbosity_, bool isBacktrack_, double rho_,
                    double c1_, double ctol_, bool isAugmentedLagrangian_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_), ctol(ctol_),
        isAugmentedLagrangian(isAugmentedLagrangian_) {
    if (verbosity > 1) {
      std::cout << "Running Conjugate Gradient propagator ..." << std::endl;
    }
    // check the validity of parameter
    checkParameters("conjugate gradient");
  }
  bool integrate();

  /**
   * @brief Conjugate Gradient stepper
   */
  void step();

  /**
   * @brief Conjugate Gradient status computation and thresholding
   */
  void status();
};

// ==========================================================
// =============      FeedForward Sweep         =============
// ==========================================================
/**
 * @brief Conjugate Gradient based feedforward parameter sweep
 * @param H__, vector of sweep parameters (H0)
 * @param VP__, vector of sweep parameters (V/cam)
 * @return
 */
class DLL_PUBLIC FeedForwardSweep : public ConjugateGradient {
public:
  std::vector<double> H_;
  std::vector<double> VP_;
  FeedForwardSweep(System &f_, double dt_, bool isAdaptiveStep_,
                   double total_time_, double tSave_, double tolerance_,
                   std::string outputDir_, std::string trajFileName_,
                   size_t verbosity_, bool isBacktrack_, double rho_,
                   double c1_, double ctol_, bool isAugmentedLagrangian_,
                   std::vector<double> H__, std::vector<double> VP__)
      : ConjugateGradient(f_, dt_, isAdaptiveStep_, total_time_, tSave_,
                          tolerance_, outputDir_, trajFileName_, verbosity_,
                          isBacktrack_, rho_, c1_, ctol_,
                          isAugmentedLagrangian_),
        H_(H__), VP_(VP__) {}
  void sweep();
};

} // namespace mem3dg
