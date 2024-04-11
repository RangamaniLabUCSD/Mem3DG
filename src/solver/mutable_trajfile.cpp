/*
 * Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
 *
 * Copyright 2020- The Mem3DG Authors
 * and the project initiators Cuncheng Zhu, Christopher T. Lee, and
 * Padmini Rangamani.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Please help us support Mem3DG development by citing the research
 * papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
 * for more information.
 */

/**
 * @file  mutable_trajfile.cpp
 * @brief Netcdf trajectory output support
 *
 */

#ifdef MEM3DG_WITH_NETCDF

#include <cassert>
#include <iostream>
#include <netcdf>
#include <vector>

#include "mem3dg/solver/mutable_trajfile.h"

namespace mem3dg {
namespace solver {
namespace nc = ::netCDF;
using NcException = nc::exceptions::NcException;
using NcFile = nc::NcFile;

bool MutableTrajFile::check_metadata() {
  // validate data
  std::string tmp;
  fd->getAtt(CONVENTIONS_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VALUE)
    mem3dg_runtime_error("NetCDF convention mismatch. This file does "
                         "not appear to be a valid Mem3DG trajectory.");

  fd->getAtt(CONVENTIONS_VERSION_NAME).getValues(tmp);
  if (tmp != CONVENTIONS_VERSION_VALUE)
    mem3dg_runtime_error(
        "Trajectory version mismatch. This file was generated with a "
        "different convention version.");

  return true;
}
} // namespace solver
} // namespace mem3dg

#endif
