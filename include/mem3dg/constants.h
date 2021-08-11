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

namespace mem3dg {
namespace constants {
const double kBoltzmann =
    1.380649e-8; // Boltzmann constant (nanonewton * micrometer / Kelvin)
const double PI = 3.14159265358979323846;
const double i = 1.0;
const double N = 6.02214076e5; // Avogadro constant (/atto-mol)
const double R =
    kBoltzmann *
    N; // ideal gas constant (nanonewton * micrometer / Kelvin / atto-mol)
} // namespace constants
} // namespace mem3dg
