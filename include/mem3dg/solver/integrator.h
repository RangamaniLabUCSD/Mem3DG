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
#include <geometrycentral/surface/rich_surface_mesh_data.h>
#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
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
namespace integration {

// ==========================================================
// =============      Time Integration         ==============
// ==========================================================
DLL_PUBLIC void stormerVerlet(System &f, double dt, double total_time,
                              double tolerance);
/**
 * @brief Velocity Verlet time Integration
 * @param f, force object
 * @param dt, time step
 * @param total_time, total simulation time
 * @param tSave, period of saving output data
 * @param tolerance, tolerance for termination
 * @param verbosity,   0: no output file 1: Traj.nc 2: All above + console
 * printing 3: All above + txt + .ply
 * @param isAdaptiveStep, option to scale time step according to mesh size
 * @param outputDir, path to the output directory
 * @return
 */
DLL_PUBLIC void velocityVerlet(System &f, double dt, double total_time,
                               double tSave, double tolerance,
                               const size_t verbosity,
                               const bool isAdaptiveStep,
                               std::string outputDir);
/**
 * @brief Euler (gradient descent) time Integration
 * @param f, force object
 * @param dt, time step
 * @param total_time, total simulation time
 * @param tSave, period of saving output data
 * @param tolerance, tolerance for termination
 * @param verbosity,   0: no output file 1: Traj.nc 2: All above + console
 * printing 3: All above + txt + .ply
 * @param isAdaptiveStep, option to scale time step according to mesh size
 * @param outputDir, path to the output directory
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @return
 */
DLL_PUBLIC bool euler(System &f, double dt, double total_time, double tSave,
                      double tolerance, const size_t verbosity,
                      std::string outputDir, const bool isBacktrack,
                      const double rho, const double c1,
                      const bool isAdaptiveStep);

/**
 * @brief Conjugate Gradient propagator
 * @param f, force object
 * @param dt, time step
 * @param total_time, total simulation time
 * @param tSave, period of saving output data
 * @param tol, tolerance for termination (force)
 * @param ctol, tolerance for termination (contraints)
 * @param verbosity,   0: no output file 1: Traj.nc 2: All above + console
 * printing 3: All above + txt + .ply
 * @param isAdaptiveStep, option to scale time step according to mesh size
 * @param outputDir, path to the output directory
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @param tragFileName, file name for output trajectory
 * @return
 */
DLL_PUBLIC bool conjugateGradient(System &f, double dt, double total_time,
                                  double tSave, double tol, double ctol,
                                  const size_t verbosity, std::string outputDir,
                                  const bool isBacktrack, const double rho,
                                  const double c1,
                                  const bool isAugmentedLagrangian,
                                  const bool isAdaptiveStep,
                                  const std::string trajFileName);

/**
 * @brief Conjugate Gradient based feedforward parameter sweep
 * @param f, force object
 * @param dt, time step
 * @param maxtime, max simulation time for each set of parameters
 * @param tSave, period of saving output data
 * @param tol, tolerance for termination (force)
 * @param ctol, tolerance for termination (contraints)
 * @param isAdaptiveStep, option to scale time step according to mesh size
 * @param outputDir, path to the output directory
 * @param isBacktrack, option to use backtracking line search algorithm
 * @param rho, backtracking coefficient
 * @param c1, Wolfe condition parameter
 * @param isAugmentedLagrangian, option to use Augmented Lagrangian method
 * @return
 */
DLL_PUBLIC void
feedForwardSweep(System &f, std::vector<double> H_, std::vector<double> V_,
                 double dt, double maxTime, double tSave, double tol,
                 double ctol, std::string outputDir, const bool isBacktrack,
                 const double rho, const double c1,
                 const bool isAugmentedLagrangian, const bool isAdaptiveStep);

// ==========================================================
// =================     Output Data         ================
// ==========================================================
/**
 * @brief Save data to richData
 * @param f, force object
 * @param physcialPressure, physical pressre eigen matrix
 * @param verbosity, verbosity setting
 * @return
 */
void saveRichData(
    System &f, const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t verbosity);

/**
 * @brief Save data to netcdf traj file
 * @param f, force object
 * @param frame, frame index of netcdf traj file
 * @param time, simulation time
 * @param fd, netcdf trajFile object
 * @param physcialPressure, physical pressre eigen matrix
 * @param energy, components of energy - totalE, BE, sE, pE, kE, cE, lE,
 * exE
 * @param verbosity, verbosity setting
 * @return
 */
void saveNetcdfData(
    const System &f, size_t &frame, const double &time, TrajFile &fd,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
    const size_t &verbosity);

/**
 * @brief Stomer Verlet time Integration
 * @param f, force object
 * @param dt, time step
 * @param total_time, total simulation time
 * @param tolerance, tolerance for termination
 * @return
 */

DLL_PUBLIC
void getParameterLog(System &f, double dt, double total_time, double tolerance,
                     double tSave, std::string inputMesh,
                     std::string outputDir);

DLL_PUBLIC void
getStatusLog(std::string nameOfFile, System &f, double dt, double time,
             std::size_t frame, double areaError, double volumeError,
             double bendingError, double faceError, double bendingEnergy,
             double surfaceEnergy, double pressureEnergy, double kineticEnergy,
             double chemicalEnergy, double lineEnergy, double totalEnergy,
             double L2ErrorNorm, bool isTuftedLaplacian, bool isProtein,
             bool isVertexShift, std::string inputMesh);

void getEnergyLog(double time, double bendingEnergy, double surfaceEnergy,
                  double pressureEnergy, double kineticEnergy,
                  double chemicalEnergy, double totalEnergy,
                  std::string outputDir);

// ==========================================================
// =============     Helper functions          ==============
// ==========================================================
/**
 * @brief Thresholding when adopting reduced volume parametrization
 * @param f, reference to the system object
 * @param EXIT, reference to the exit flag
 * @param isAugmentedLagrangian, whether using augmented lagrangian method
 * @param dArea, normalized area difference
 * @param dVolume, normalized volume difference
 * @param ctol, exit criterion for constraint
 * @param tol, exit criterion for gradient
 * @param increment, increment coefficient of penalty when using incremental
 * penalty method
 * @return
 */
void reducedVolumeThreshold(System &f, bool &EXIT,
                            const bool isAugmentedLagrangian,
                            const double dArea, const double dVolume,
                            const double ctol, const double tol,
                            double increment);
/**
 * @brief Thresholding when adopting ambient pressure constraint
 * @param f, reference to the system object
 * @param EXIT, reference to the exit flag
 * @param isAugmentedLagrangian, whether using augmented lagrangian method
 * @param dArea, normalized area difference
 * @param ctol, exit criterion for constraint
 * @param tol, exit criterion for gradient
 * @param increment, increment coefficient of penalty when using incremental
 * penalty method
 * @return
 */
void pressureConstraintThreshold(System &f, bool &EXIT,
                                 const bool isAugmentedLagrangian,
                                 const double dArea, const double ctol,
                                 const double tol, double increment);

/**
 * @brief Backtracking algorithm that dynamically adjust step size based on
 * energy evaluation
 * @param f, force object
 * @param dt, initial step size
 * @param rho, discount factor
 * @param c1, constant for Wolfe condtion, between 0 to 1, usually ~ 1e-4
 * @param time, simulation time
 * @param EXIT, exit flag for integration loop
 * @param SUCCESS, test flag on simulation sucess
 * @param verbosity, verbosity setting
 * @param potentialEnergy_pre, previous energy evaluation
 * @param force, gradient of the energy
 * @param direction, direction, most likely some function of gradient
 * @return
 */
void backtrack(System &f, const double dt, double rho, double c1, double &time,
               bool &EXIT, bool &SUCCESS, const size_t verbosity,
               const double potentialEnergy_pre,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &force,
               const Eigen::Matrix<double, Eigen::Dynamic, 3> &direction);

/**
 * @brief Summerize forces into 3 categories: physcialPressure, DPDPressure and
 * regularizationForce. Note that the forces has been removed rigid body mode
 * and masked for integration
 *
 * @param f
 * @param physicalPressure
 * @param DPDPressure
 * @param regularizationForce
 * @return
 */
void getForces(System &f,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &physicalPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &DPDPressure,
               Eigen::Matrix<double, Eigen::Dynamic, 3> &regularizationForce);

/**
 * @brief Check parameters for time integration
 *
 * @param integrator, name of time integration
 * @param f, system object
 * @param h, time step
 * @return
 */
void checkParameters(std::string integrator, System &f, double h);


} // namespace integration
} // namespace mem3dg
