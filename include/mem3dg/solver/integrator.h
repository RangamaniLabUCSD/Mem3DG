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
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

#include "mem3dg/solver/meshops.h"
#include "system.h"
#include <csignal>
#include <stdexcept>
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
  /// Energy of the last time step
  Energy previousE;
  /// time step
  double dt;
  // total simulation time
  double total_time;
  /// period of saving output data
  double tSave;
  /// last time saving the data
  double lastSave;
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
  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalForceVec;
  Eigen::Matrix<double, Eigen::Dynamic, 1> physicalForce;
  /// numerical dissipative particle dynamics force to the system
  Eigen::Matrix<double, Eigen::Dynamic, 1> DPDForce;

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
  /// initial maximum force
  double maxForce;

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
      : f(f_), previousE(f_.E), dt(dt_), isAdaptiveStep(isAdaptiveStep_),
        total_time(total_time_), tSave(tSave_), tol(tolerance_),
        verbosity(verbosity_), outputDir(outputDir_),
        trajFileName(trajFileName_), init_time(f.time), SUCCESS(true),
        EXIT(false), frame(0) {

    // Initialize the timestep-meshsize ratio
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();

    // Initialize the initial maxForce
    getForces();
    maxForce = f.O.isShapeVariation
                   ? physicalForce.cwiseAbs().maxCoeff()
                   : f.F.chemicalPotential.raw().cwiseAbs().maxCoeff();

    // Initialize geometry constraints
    dArea = 1e10;
    dVP = 1e10;

    // Initialize system summarized forces
    regularizationForce.resize(f.mesh->nVertices(), 3);
    physicalForce.resize(f.mesh->nVertices(), 1);
    DPDForce.resize(f.mesh->nVertices(), 1);
    regularizationForce.setZero();
    physicalForce.setZero();
    DPDForce.setZero();
  }

  // ==========================================================
  // =================   Template functions    ================
  // ==========================================================
  // virtual bool integrate() { return true; };
  // virtual void march(){};
  // virtual void status(){};
  // virtual void checkParameters(){};

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
   * @brief Save data to richData and output ply file
   */
  void saveRichData(std::string plyName);

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Save data to netcdf traj file
   * @param frame, frame index of netcdf traj file
   * @param fd, netcdf trajFile object
   */
  void saveNetcdfData();

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
                                const EigenVectorX3D &initial_pos,
                                const EigenVectorX1D &init_proteinDensity,
                                bool runAll = false);
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
                   outputDir_, trajFileName_, verbosity_) {
    throw std::runtime_error(
        "StomerVerlet is currently not tested and maintained!");
  }

  bool integrate();
};

// ==========================================================
// =============        Velocity Verlet        +=============
// ==========================================================
/**
 * @brief Velocity Verlet time Integration
 */
class DLL_PUBLIC VelocityVerlet : public Integrator {
public:
  // total pressure
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure;
  // total pressure of new iteration
  Eigen::Matrix<double, Eigen::Dynamic, 3> newTotalPressure;
  // total energy of the system
  double totalEnergy;

  VelocityVerlet(System &f_, double dt_, bool isAdaptiveStep_,
                 double total_time_, double tSave_, double tolerance_,
                 std::string outputDir_, std::string trajFileName_,
                 size_t verbosity_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_) {

    // print to console
    if (verbosity > 1) {
      std::cout << "Running Velocity Verlet integrator ..." << std::endl;
    }

    // check the validity of parameter
    checkParameters();

    totalPressure.resize(f.mesh->nVertices(), 3);
    newTotalPressure.resize(f.mesh->nVertices(), 3);
    totalPressure.setZero();
    newTotalPressure.setZero();

    f.computeFreeEnergy();
    totalEnergy = f.E.totalE;
  }

  /**
   * @brief velocity Verlet driver function
   */
  bool integrate();
  /**
   * @brief velocity Verlet marcher
   */
  void march();

  /**
   * @brief Velocity Verlet status computation and thresholding
   */
  void status();

  /**
   * @brief Check parameters for time integration
   */
  void checkParameters();

  /**
   * @brief step for n iterations
   */
  void step(size_t n) {
    for (size_t i = 0; i < n; i++) {
      status();
      march();
    }
  }
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

    // print to console
    if (verbosity > 1) {
      std::cout << "Running Forward Euler (steepest descent) propagator ..."
                << std::endl;
    }
    // check the validity of parameter
    checkParameters();
  }

  /**
   * @brief Forward Euler driver function
   */
  bool integrate();

  /**
   * @brief Forward Euler stepper
   */
  void march();

  /**
   * @brief Forward Euler status computation and thresholding
   */
  void status();

  /**
   * @brief Check parameters for time integration
   */
  void checkParameters();

  /**
   * @brief step for n iterations
   */
  void step(size_t n) {
    for (size_t i = 0; i < n; i++) {
      status();
      march();
    }
  }
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

  double currentNormSq;
  double pastNormSq;

  size_t countCG = 0;
  size_t restartNum;

  // size_t countPM = 0;

  ConjugateGradient(System &f_, double dt_, bool isAdaptiveStep_,
                    double total_time_, double tSave_, double tolerance_,
                    std::string outputDir_, std::string trajFileName_,
                    size_t verbosity_, bool isBacktrack_, double rho_,
                    double c1_, double ctol_, bool isAugmentedLagrangian_,
                    size_t restartNum_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_), ctol(ctol_),
        isAugmentedLagrangian(isAugmentedLagrangian_), restartNum(restartNum_) {

    // print to console
    if (verbosity > 1) {
      std::cout << "Running Conjugate Gradient propagator ..." << std::endl;
    }

    // check the validity of parameter
    checkParameters();
  }

  /**
   * @brief Conjugate Gradient driver function
   */
  bool integrate();

  /**
   * @brief Conjugate Gradient stepper
   */
  void march();

  /**
   * @brief Conjugate Gradient status computation and thresholding
   */
  void status();

  /**
   * @brief Check parameters for time integration
   */
  void checkParameters();

  /**
   * @brief step for n iterations
   */
  void step(size_t n) {
    for (size_t i = 0; i < n; i++) {
      status();
      march();
    }
  }
};

// ==========================================================
// =============              BFGS              =============
// ==========================================================
/**
 * @brief BFGS optimizer
 * @param ctol, tolerance for termination (contraints)
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @return Success, if simulation is sucessful
 */
class DLL_PUBLIC BFGS : public Integrator {
public:
  const bool isBacktrack;
  const double rho;
  const double c1;
  double ctol;
  const bool isAugmentedLagrangian;
  bool ifRestart = false;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s_protein;

  BFGS(System &f_, double dt_, bool isAdaptiveStep_, double total_time_,
       double tSave_, double tolerance_, std::string outputDir_,
       std::string trajFileName_, size_t verbosity_, bool isBacktrack_,
       double rho_, double c1_, double ctol_, bool isAugmentedLagrangian_)
      : Integrator(f_, dt_, isAdaptiveStep_, total_time_, tSave_, tolerance_,
                   outputDir_, trajFileName_, verbosity_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_), ctol(ctol_),
        isAugmentedLagrangian(isAugmentedLagrangian_) {

    throw std::runtime_error("BFGS is currently not tested and maintained!");

    // print to console
    if (verbosity > 1) {
      std::cout << "Running BFGS propagator ..." << std::endl;
    }

    hess_inv.resize(f.mesh->nVertices() * 3, f.mesh->nVertices() * 3);
    hess_inv.setIdentity();
    pastPhysicalForce.resize(f.mesh->nVertices() * 3, 1);
    pastPhysicalForce.setZero();
    s.resize(f.mesh->nVertices() * 3, 1);
    s.setZero();

    hess_inv_protein.resize(f.mesh->nVertices(), f.mesh->nVertices());
    hess_inv_protein.setIdentity();
    pastPhysicalForce_protein.resize(f.mesh->nVertices(), 1);
    pastPhysicalForce_protein.setZero();
    s_protein.resize(f.mesh->nVertices(), 1);
    s_protein.setZero();

    // check the validity of parameter
    checkParameters();
  }

  /**
   * @brief BFGS function
   */
  bool integrate();

  /**
   * @brief BFGS stepper
   */
  void march();

  /**
   * @brief BFGS status computation and thresholding
   */
  void status();

  /**
   * @brief Check parameters for time integration
   */
  void checkParameters();

  /**
   * @brief step for n iterations
   */
  void step(size_t n) {
    for (size_t i = 0; i < n; i++) {
      status();
      march();
    }
  }
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
                   size_t restartNum_, std::vector<double> H__,
                   std::vector<double> VP__)
      : ConjugateGradient(f_, dt_, isAdaptiveStep_, total_time_, tSave_,
                          tolerance_, outputDir_, trajFileName_, verbosity_,
                          isBacktrack_, rho_, c1_, ctol_,
                          isAugmentedLagrangian_, restartNum_),
        H_(H__), VP_(VP__) {
    throw std::runtime_error(
        "FeedForwardSweep is currently not tested and maintained!");
  }
  void sweep();
};

} // namespace mem3dg
