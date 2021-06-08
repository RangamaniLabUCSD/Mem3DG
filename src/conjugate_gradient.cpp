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
#include "mem3dg/solver/integrator.h"
#include "mem3dg/solver/meshops.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
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
  createNetcdfFile();
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
      f.processMesh();
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
  if (f.P.gamma != 0 || f.P.temp != 0) {
    throw std::runtime_error("DPD has to be turned off for CG integration!");
  }
  if (f.P.Bc != 1 && f.P.Bc != 0) {
    throw std::runtime_error(
        "Binding constant should be set to 1 for optimization!");
  }
  if (isBacktrack) {
    if (rho >= 1 || rho <= 0 || c1 >= 1 || c1 <= 0) {
      throw std::runtime_error("To backtrack, 0<rho<1 and 0<c1<1!");
    }
  }
  if (restartNum < 1){
    throw std::runtime_error("restartNum > 0!");
  }
  // if (f.O.isVertexShift) {
  //   throw std::runtime_error(
  //       "Vertex shift is not supported for CG integration!");
  // }
}

void ConjugateGradient::status() {
  auto physicalForce = f.F.toMatrix(f.F.mechanicalForce);

  // compute summerized forces
  getForces();

  // compute the area contraint error
  dArea = abs(f.surfaceArea / f.refSurfaceArea - 1);
  if (f.O.isReducedVolume) {
    dVP = abs(f.volume / f.refVolume / f.P.Vt - 1);
    reducedVolumeThreshold(EXIT, isAugmentedLagrangian, dArea, dVP, ctol, 1.3);
  } else {
    dVP = abs(1 / f.volume / f.P.cam - 1);
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
  auto vel_e = f.F.toMatrix(f.vel);
  auto vel_protein_e = f.F.toMatrix(f.vel_protein);
  auto pos_e = f.F.toMatrix(f.vpg->inputVertexPositions);
  auto physicalForceVec = f.F.toMatrix(f.F.mechanicalForceVec);
  auto physicalForce = f.F.toMatrix(f.F.mechanicalForce);
  // typedef gc::EigenVectorMap_T<double, 3,
  // Eigen::RowMajor>(*EigenMap3)(gcs::VertexData<gc::Vector3>); EigenMap3
  // Map3 = gc::EigenMap<double, 3>;

  // determine conjugate gradient direction, restart after nVertices() cycles
  if (countCG % restartNum == 0) {
    pastNormSq =
        (f.O.isShapeVariation ? physicalForceVec.squaredNorm() : 0) +
        (f.O.isProteinVariation ? f.F.chemicalPotential.raw().squaredNorm()
                                : 0);
    vel_e = physicalForceVec;
    vel_protein_e = f.F.chemicalPotential.raw();
    countCG = 1;
  } else {
    currentNormSq =
        (f.O.isShapeVariation ? physicalForceVec.squaredNorm() : 0) +
        (f.O.isProteinVariation ? f.F.chemicalPotential.raw().squaredNorm()
                                : 0);
    vel_e *= currentNormSq / pastNormSq;
    vel_e += physicalForceVec;
    vel_protein_e *= currentNormSq / pastNormSq;
    vel_protein_e += f.F.chemicalPotential.raw();
    pastNormSq = currentNormSq;
    countCG++;
  }

  // adjust time step if adopt adaptive time step based on mesh size
  if (isAdaptiveStep) {
    double minMeshLength = f.vpg->edgeLengths.raw().minCoeff();
    dt = dt_size2_ratio * maxForce * minMeshLength * minMeshLength /
         (f.O.isShapeVariation
              ? physicalForce.cwiseAbs().maxCoeff()
              : f.F.chemicalPotential.raw().cwiseAbs().maxCoeff());
  }

  // time stepping on vertex position
  previousE = f.E;
  if (isBacktrack) {
    backtrack(f.E.potE, vel_e, vel_protein_e, rho, c1);
  } else {
    pos_e += vel_e * dt;
    f.proteinDensity.raw() += f.P.Bc * vel_protein_e * dt;
    f.time += dt;
  }

  // regularization
  if ((f.P.Kse != 0) || (f.P.Ksl != 0) || (f.P.Kst != 0)) {
    f.computeRegularizationForce();
    f.vpg->inputVertexPositions.raw() += f.F.regularizationForce.raw();
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
  const double KV = f.P.Kv, KSG = f.P.Ksg, init_time = 0.0;
  const size_t verbosity = 2;

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
        f.P.lambdaSG = 0;
        f.P.lambdaV = 0;
      } else {
        f.P.Kv = KV;
        f.P.Ksg = KSG;
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
      f.P.H0c = H;
      (f.O.isReducedVolume ? f.P.Vt : f.P.cam) = VP;
      std::cout << "\nH0c: " << f.P.H0c << std::endl;
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

} // namespace mem3dg
