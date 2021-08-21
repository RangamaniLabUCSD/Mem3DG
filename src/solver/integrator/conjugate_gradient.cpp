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

#include <Eigen/Core>
#include <iostream>
#include <math.h>
#include <pcg_random.hpp>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#include <geometrycentral/utilities/vector3.h>
#include <stdexcept>

#include "Eigen/src/Core/util/Constants.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "mem3dg/meshops.h"
#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
namespace gc = ::geometrycentral;

bool ConjugateGradient::integrate() {

  signal(SIGINT, signalHandler);

#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // initialize netcdf traj file
#ifdef MEM3DG_WITH_NETCDF
  // createNetcdfFile();
  createMutableNetcdfFile();
  // print to console
  std::cout << "Initialized NetCDF file at " << outputDir + "/" + trajFileName
            << std::endl;
#endif

  // time integration loop
  for (;;) {

    // Evaluate and threhold status data
    status();

    // Save files every tSave period and print some info
    if (f.time - lastSave >= tSave || f.time == init_time || EXIT) {
      lastSave = f.time;
      saveData();
    }

    // Process mesh every tProcessMesh period
    if (f.time - lastProcessMesh > tProcessMesh) {
      lastProcessMesh = f.time;
      f.mutateMesh();
      f.updateVertexPositions(false);
    }

    // update geodesics every tUpdateGeodesics period
    if (f.time - lastUpdateGeodesics > tUpdateGeodesics) {
      lastUpdateGeodesics = f.time;
      f.updateVertexPositions(true);
    }

    // break loop if EXIT flag is on
    if (EXIT) {
      break;
    }

    // step forward
    if (f.time == lastProcessMesh || f.time == lastUpdateGeodesics) {
      f.time += 1e-10 * dt;
      countCG = 0;
    } else {
      march();
    }
  }

  // return if optimization is sucessful
  if (!SUCCESS) {
    if (tol == 0) {
      markFileName("_most");
    } else {
      markFileName("_failed");
    }
  }

  // stop the timer and report time spent
#ifdef __linux__
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal integration time: " << duration << " seconds"
              << std::endl;
  }
#endif

  return SUCCESS;
}

void ConjugateGradient::checkParameters() {
  if (f.parameters.dpd.gamma != 0) {
    mem3dg_runtime_error("DPD has to be turned off for CG integration!");
  }
  if (f.parameters.proteinMobility != 1 && f.parameters.proteinMobility != 0) {
    mem3dg_runtime_error("Protein mobility constant should "
                         "be set to 1 for optimization!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      mem3dg_runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
  if (restartNum < 1) {
    mem3dg_runtime_error("restartNum > 0!");
  }
  // if (f.O.isVertexShift) {
  //   mem3dg_runtime_error(
  //       "Vertex shift is not supported for CG integration!");
  // }
}

void ConjugateGradient::status() {
  auto physicalForce = toMatrix(f.forces.mechanicalForce);

  // compute summerized forces
  getForces();

  // compute the area contraint error
  dArea = abs(f.surfaceArea / f.parameters.tension.At - 1);
  if (f.parameters.osmotic.isPreferredVolume) {
    dVP = abs(f.volume / f.parameters.osmotic.Vt - 1);
    reducedVolumeThreshold(EXIT, isAugmentedLagrangian, dArea, dVP, ctol, 1.3);
  } else {
    dVP =
        abs(f.parameters.osmotic.n / f.volume / f.parameters.osmotic.cam - 1.0);
    pressureConstraintThreshold(EXIT, isAugmentedLagrangian, dArea, ctol, 1.3);
  }

  // exit if reached time
  if (f.time > total_time) {
    std::cout << "\nReached time." << std::endl;
    EXIT = true;
    SUCCESS = false;
  }

  // compute the free energy of the system
  f.computeFreeEnergy();

  // backtracking for error
  finitenessErrorBacktrack();
}

void ConjugateGradient::march() {
  // map the raw eigen datatype for computation
  auto vel_e = toMatrix(f.velocity);
  auto vel_protein_e = toMatrix(f.proteinVelocity);
  auto pos_e = toMatrix(f.vpg->inputVertexPositions);
  auto physicalForceVec = toMatrix(f.forces.mechanicalForceVec);
  auto physicalForce = toMatrix(f.forces.mechanicalForce);
  // typedef gc::EigenVectorMap_T<double, 3,
  // Eigen::RowMajor>(*EigenMap3)(gcs::VertexData<gc::Vector3>); EigenMap3
  // Map3 = gc::EigenMap<double, 3>;

  // determine conjugate gradient direction, restart after nVertices() cycles
  if (countCG % restartNum == 0) {
    pastNormSq = (f.parameters.variation.isShapeVariation
                      ? physicalForceVec.squaredNorm()
                      : 0) +
                 (f.parameters.variation.isProteinVariation
                      ? f.forces.chemicalPotential.raw().squaredNorm()
                      : 0);
    vel_e = physicalForceVec;
    vel_protein_e =
        f.parameters.proteinMobility * f.forces.chemicalPotential.raw();
    countCG = 1;
  } else {
    currentNormSq = (f.parameters.variation.isShapeVariation
                         ? physicalForceVec.squaredNorm()
                         : 0) +
                    (f.parameters.variation.isProteinVariation
                         ? f.forces.chemicalPotential.raw().squaredNorm()
                         : 0);
    vel_e *= currentNormSq / pastNormSq;
    vel_e += physicalForceVec;
    vel_protein_e *= currentNormSq / pastNormSq;
    vel_protein_e +=
        f.parameters.proteinMobility * f.forces.chemicalPotential.raw();
    pastNormSq = currentNormSq;
    countCG++;
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
    dt = dt_size2_ratio * maxForce * minMeshLength * minMeshLength /
         (f.parameters.variation.isShapeVariation
              ? physicalForce.cwiseAbs().maxCoeff()
              : vel_protein_e.cwiseAbs().maxCoeff());
  }

  // time stepping on vertex position
  previousE = f.energy;
  if (isBacktrack) {
    backtrack(f.energy.potE, vel_e, vel_protein_e, rho, c1);
  } else {
    pos_e += vel_e * dt;
    f.proteinDensity.raw() += vel_protein_e * dt;
    f.time += dt;
  }

  // regularization
  if (f.meshProcessor.isMeshRegularize) {
    f.computeRegularizationForce();
    f.vpg->inputVertexPositions.raw() += f.forces.regularizationForce.raw();
  }

  // recompute cached values
  f.updateVertexPositions(false);
}

void FeedForwardSweep::sweep() {
#ifdef __linux__
  // start the timer
  struct timeval start;
  gettimeofday(&start, NULL);
#endif

  // initialize variables
  const double KV = f.parameters.osmotic.Kv, KSG = f.parameters.tension.Ksg,
               init_time = 0.0;
  const std::size_t verbosity = 2;

  // initialize variables used if adopting adaptive time step based on mesh
  // size
  double dt_size2_ratio;
  if (isAdaptiveStep) {
    dt_size2_ratio = dt / f.vpg->edgeLengths.raw().minCoeff() /
                     f.vpg->edgeLengths.raw().minCoeff();
  }

  // parameter sweep
  for (double H : H_) {
    for (double VP : VP_) {
      // reset parameters
      if (isAugmentedLagrangian) {
        f.parameters.tension.lambdaSG = 0;
        f.parameters.osmotic.lambdaV = 0;
      } else {
        f.parameters.osmotic.Kv = KV;
        f.parameters.tension.Ksg = KSG;
      }

      // adjust time step if adopt adaptive time step based on mesh size
      if (isAdaptiveStep) {
        double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
        dt = dt_size2_ratio * minMeshLength * minMeshLength;
      }

      // initialize trajectory file name
      char buffer[50];
      sprintf(buffer, "traj_H_%d_VP_%d.nc", int(H * 100), int(VP * 100));
      trajFileName = buffer;

      // update sweeping paraemters
      f.parameters.bending.H0c = H;
      (f.parameters.osmotic.isPreferredVolume ? f.parameters.osmotic.Vt
                                              : f.parameters.osmotic.cam) = VP;
      std::cout << "\nH0c: " << f.parameters.bending.H0c << std::endl;
      std::cout << "VP: " << VP << std::endl;

      // recalculate cached values
      f.updateVertexPositions();

      // rerun CG optimization
      integrate();
    }
    // reverse the order of inner loop to ensure phase space closeness
    std::reverse(VP_.begin(), VP_.end());
  }

#ifdef __linux__
  // stop the timer and report time spent
  double duration = getDuration(start);
  if (verbosity > 0) {
    std::cout << "\nTotal parameter sweep time: " << duration << " seconds"
              << std::endl;
  }
#endif
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
