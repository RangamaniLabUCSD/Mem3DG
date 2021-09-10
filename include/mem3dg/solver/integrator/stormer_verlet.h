// Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) 2021:
//     Laboratory for Computational Cellular Mechanobiology
//     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
//     Christopher T. Lee (ctlee@ucsd.edu)
//     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
//     Padmini Rangamani (prangamani@eng.ucsd.edu)
//

#pragma once

#include "mem3dg/solver/integrator/integrator.h"
#include "mem3dg/solver/system.h"

namespace mem3dg {
namespace solver {
namespace integrator {
/**
 * @brief Stomer Verlet time Integration
 */
class DLL_PUBLIC StormerVerlet : public Integrator {
private:
  gcs::VertexData<gc::Vector3> pastPositions;

public:
  StormerVerlet(System &system_, double characteristicTimeStep_,
                double totalTime_, double savePeriod_, double tolerance_,
                std::string outputDirectory_)
      : Integrator(system_, characteristicTimeStep_, totalTime_, savePeriod_,
                   tolerance_, outputDirectory_),
        pastPositions(*system.mesh, {0, 0, 0}) {
    mem3dg_runtime_error(
        "StomerVerlet is currently not tested and maintained!");
  }

  bool integrate() override;
  void march() override{};
  void status() override{};
  void checkParameters() override{};
};
} // namespace integrator
} // namespace solver
} // namespace mem3dg
