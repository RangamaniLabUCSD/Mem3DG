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
class DLL_PUBLIC Integrator {
public:
  System &f;
  double dt;
  double total_time;
  double tSave;
  size_t verbosity;
  double tol;
  std::string outputDir;
  bool isAdaptiveStep;
  std::string trajFileName;
  Eigen::Matrix<double, Eigen::Dynamic, 3> regularizationForce;
  Eigen::Matrix<double, Eigen::Dynamic, 3> physicalPressure;
  Eigen::Matrix<double, Eigen::Dynamic, 3> DPDPressure;

  /**
   * @brief Construct a new integrator object
   * @param f, force object
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
             std::string trajFileName_ = "/traj.nc", size_t verbosity_ = 3)
      : f(f_), dt(dt_), isAdaptiveStep(isAdaptiveStep_),
        total_time(total_time_), tSave(tSave_), tol(tolerance_),
        verbosity(verbosity_), outputDir(outputDir_),
        trajFileName(trajFileName_) {
    std::cout << "Solving the system and saving to " << outputDir << std::endl;
  }

  // ==========================================================
  // =============      Time Integration         ==============
  // ==========================================================
  /**
   * @brief Stomer Verlet time Integration
   * @return
   */
  void stormerVerlet();

  /**
   * @brief Velocity Verlet time Integration
   * @return
   */
  void velocityVerlet();

  /**
   * @brief Euler (gradient descent) time Integration
   * @param isBacktrack, option to use backtracking line search algorithm
   * @param rho, backtracking coefficient
   * @param c1, Wolfe condition parameter
   * @return
   */
  bool euler(const bool isBacktrack, const double rho, const double c1);

  /**
   * @brief Conjugate Gradient propagator
   * @param ctol, tolerance for termination (contraints)
   * @param isBacktrack, option to use backtracking line search algorithm
   * @param rho, backtracking coefficient
   * @param c1, Wolfe condition parameter
   * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
   * @return
   */
  bool conjugateGradient(double ctol, const bool isBacktrack, const double rho,
                         const double c1, const bool isAugmentedLagrangian);

  /**
   * @brief Conjugate Gradient based feedforward parameter sweep
   * @param ctol, tolerance for termination (contraints)
   * @param isBacktrack, option to use backtracking line search algorithm
   * @param rho, backtracking coefficient
   * @param c1, Wolfe condition parameter
   * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
   * @return
   */
  void feedForwardSweep(std::vector<double> H_, std::vector<double> VP_,
                        double ctol, const bool isBacktrack, const double rho,
                        const double c1, const bool isAugmentedLagrangian);

  // ==========================================================
  // =================     Output Data         ================
  // ==========================================================
  /**
   * @brief Save data to richData
   * @return
   */
  void saveRichData();

#ifdef MEM3DG_WITH_NETCDF
  /**
   * @brief Save data to netcdf traj file
   * @param frame, frame index of netcdf traj file
   * @param fd, netcdf trajFile object
   * @return
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

} // namespace mem3dg
