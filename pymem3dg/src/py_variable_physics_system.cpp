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

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mem3dg/mem3dg"

namespace py = pybind11;

namespace mem3dg {
namespace solver {

void init_variable_physics_system(py::module &mod) {
  py::class_<VariablePhysicsSystem> variable_physics_system(
      mod, "Variable_Physics_System",
      R"delim(
            System with variable physics
        )delim");

  variable_physics_system.def(py::init<std::list<Physic *>>(), "Constructor");

  class PyPhysic : public Physic {
  public:
    /* Inherit the constructors */
    using Physic::Physic;

    void compute(VariablePhysicsSystem& S) override {
      PYBIND11_OVERRIDE_PURE(
          void,   /* Return type */
          Physic, /* Parent class */
          compute, /* Name of function in C++ (must match Python name) */
          S
      );
    }
    void initialize(VariablePhysicsSystem& S) override {
      PYBIND11_OVERRIDE_PURE(
          void,      /* Return type */
          Physic,    /* Parent class */
          initialize, /* Name of function in C++ (must match Python name) */
          S
      );
    }
  };

  py::class_<Physic, PyPhysic> physic(mod, "Physic",
                                      R"delim(
      Physic ABC
  )delim");

  py::class_<BendingPhysic, Physic> bending_physic(mod, "Bending_Physic",
                                                   R"delim(
      Bending physic
  )delim");

  bending_physic.def(py::init<>(), "Default constructor");

  py::class_<StretchingPhysic, Physic> streching_physic(mod,
                                                        "Stretching_Physic",
                                                        R"delim(
      Stretching physic
  )delim");

  streching_physic.def(py::init<>(), "Default constructor");
}

} // namespace solver
} // namespace mem3dg
