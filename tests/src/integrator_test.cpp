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

#include <iostream>

#include <gtest/gtest.h>

#include "mem3dg/constants.h"
#include "mem3dg/mem3dg"
#include <Eigen/Core>

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

class IntegratorTest : public ::testing::Test {
public:
  IntegratorTest() {
    std::tie(mesh, vpg) = mem3dg::getIcosphereMatrix(1, 3);

    /// physical parameters
    p.bending.Kbc = 8.22e-5;

    auto preferredAreaSurfaceTensionModel = [](double area) {
      double Ksg = 0.1;
      double At = 4.0 * mem3dg::constants::PI;
      double A_difference = area - At;
      double tension = Ksg * A_difference / At;
      double energy = tension * A_difference / 2;
      return std::make_tuple(tension, energy);
    };
    p.tension.form = preferredAreaSurfaceTensionModel;

    auto preferredVolumeOsmoticPressureModel = [](double volume) {
      double isPreferredVolume = true;
      double Kv = 0.01;
      double Vt = 4.0 / 3.0 * mem3dg::constants::PI * 0.7;
      double osmoticPressure = -(Kv * (volume - Vt) / Vt / Vt);

      double V_difference = volume - Vt;
      double pressureEnergy = -osmoticPressure * V_difference / 2;

      return std::make_tuple(osmoticPressure, pressureEnergy);
    };

    p.osmotic.form = preferredVolumeOsmoticPressureModel;
  }
  void SetUp() override {}
  void TearDown() override {}

  Eigen::Matrix<std::size_t, Eigen::Dynamic, 3, Eigen::RowMajor> mesh;
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vpg;

  mem3dg::solver::Parameters p;
  std::string trajName = "traj.nc";

  const double dt = 0.5, T = 50, eps = 0, tSave = 10;
  const std::string outputDir = "/tmp";
};

TEST_F(IntegratorTest, EulerIntegratorTest) {
  mem3dg::solver::Geometry geometry(mesh, vpg);
  mem3dg::solver::System f(geometry, p, 0);
  f.initialize(false);
  mem3dg::solver::integrator::Euler integrator{f,   dt,        T, tSave,
                                               eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.ifOutputTrajFile = true;
  integrator.integrate();
}

TEST_F(IntegratorTest, ConjugateGradientIntegratorTest) {
  mem3dg::solver::Geometry geometry(mesh, vpg);
  mem3dg::solver::System f(geometry, p, 0);
  f.initialize(false);
  mem3dg::solver::integrator::ConjugateGradient integrator{
      f, dt, T, tSave, eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.integrate();
}

TEST_F(IntegratorTest, VelocityVerletIntegratorTest) {
  p.damping = 0.1 / dt;
  mem3dg::solver::Geometry geometry(mesh, vpg);
  mem3dg::solver::System f(geometry, p, 0);
  f.initialize(false);
  mem3dg::solver::integrator::VelocityVerlet integrator{
      f, dt, 1, tSave, eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.integrate();
}

TEST_F(IntegratorTest, ContinuationEulerIntegratorTest) {
  mem3dg::solver::Geometry g(outputDir + "/" + trajName, 2);
  mem3dg::solver::System f(g, outputDir + "/" + trajName, 2, p);
  f.initialize(false);
  mem3dg::solver::integrator::Euler integrator{f,   dt,        T * 2, tSave,
                                               eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.integrate();
}
TEST_F(IntegratorTest, ContinuationConjugateGradientIntegratorTest) {
  mem3dg::solver::Geometry g(outputDir + "/" + trajName, 2);
  mem3dg::solver::System f(g, outputDir + "/" + trajName, 2, p);
  f.initialize(false);
  mem3dg::solver::integrator::ConjugateGradient integrator{
      f, dt, T * 2, tSave, eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.integrate();
}
TEST_F(IntegratorTest, ContinuationVelocityVerletIntegratorTest) {
  p.damping = 0.1 / dt;
  mem3dg::solver::Geometry g(outputDir + "/" + trajName, 2);
  mem3dg::solver::System f(g, outputDir + "/" + trajName, 2, p);
  f.initialize(false);
  mem3dg::solver::integrator::VelocityVerlet integrator{
      f, dt, T, tSave, eps, outputDir, 0};
  integrator.trajFileName = trajName;
  integrator.integrate();
}
