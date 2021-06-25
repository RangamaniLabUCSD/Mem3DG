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
  size_t frame;
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
  TrajFile fd;
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
      : f(f_), previousE(f_.E), dt(dt_), total_time(total_time_), tSave(tSave_),
        tol(tolerance_), tUpdateGeodesics(total_time), tProcessMesh(total_time),
        outputDir(outputDir_), init_time(f_.time), SUCCESS(true), EXIT(false),
        frame(0), lastUpdateGeodesics(f_.time), lastProcessMesh(f_.time),
        lastSave(f_.time) {

    // Initialize the timestep-meshsize ratio
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();

    // Initialize the initial maxForce
    getForces();
    maxForce = f.O.isShapeVariation
                   ? f.F.toMatrix(f.F.mechanicalForce).cwiseAbs().maxCoeff()
                   : f.F.chemicalPotential.raw().cwiseAbs().maxCoeff();

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
private:
  gcs::VertexData<gc::Vector3> pastPositions;

public:
  StormerVerlet(System &f_, double dt_, double total_time_, double tSave_,
                double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_),
        pastPositions(*f.mesh, {0, 0, 0}) {
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
private:
  // total pressure
  Eigen::Matrix<double, Eigen::Dynamic, 3> totalPressure;
  // total pressure of new iteration
  Eigen::Matrix<double, Eigen::Dynamic, 3> newTotalPressure;
  // total energy of the system
  double totalEnergy;

public:
  VelocityVerlet(System &f_, double dt_, double total_time_, double tSave_,
                 double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_) {

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
  bool isBacktrack = true;
  double rho = 0.99;
  double c1 = 0.0001;

  Euler(System &f_, double dt_, double total_time_, double tSave_,
        double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_) {

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
private:
  double currentNormSq;
  double pastNormSq;

  size_t countCG = 0;

public:
  size_t restartNum = 5;
  bool isBacktrack = true;
  double rho = 0.99;
  double c1 = 0.001;
  double ctol = 0.01;
  bool isAugmentedLagrangian = false;

  // size_t countPM = 0;

  ConjugateGradient(System &f_, double dt_, double total_time_, double tSave_,
                    double tolerance_, std::string outputDir_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_) {

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
private:
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_inv_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> pastPhysicalForce_protein;
  Eigen::Matrix<double, Eigen::Dynamic, 1> s_protein;

public:
  const bool isBacktrack;
  const double rho;
  const double c1;
  const double ctol;
  const bool isAugmentedLagrangian;
  bool ifRestart;

  BFGS(System &f_, double dt_, double total_time_, double tSave_,
       double tolerance_, std::string outputDir_, bool isAdaptiveStep_,
       std::string trajFileName_, size_t verbosity_, bool isBacktrack_,
       double rho_, double c1_, double ctol_, bool isAugmentedLagrangian_)
      : Integrator(f_, dt_, total_time_, tSave_, tolerance_, outputDir_),
        isBacktrack(isBacktrack_), rho(rho_), c1(c1_), ctol(ctol_),
        isAugmentedLagrangian(isAugmentedLagrangian_), ifRestart(false) {

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
  FeedForwardSweep(System &f_, double dt_, double total_time_, double tSave_,
                   double tolerance_, std::string outputDir_,
                   bool isAdaptiveStep_, std::string trajFileName_,
                   size_t verbosity_, bool isBacktrack_, double rho_,
                   double c1_, double ctol_, bool isAugmentedLagrangian_,
                   size_t restartNum_, std::vector<double> H__,
                   std::vector<double> VP__)
      : ConjugateGradient(f_, dt_, total_time_, tSave_, tolerance_, outputDir_),
        H_(H__), VP_(VP__) {
    throw std::runtime_error(
        "FeedForwardSweep is currently not tested and maintained!");
  }
  void sweep();
};

} // namespace mem3dg
