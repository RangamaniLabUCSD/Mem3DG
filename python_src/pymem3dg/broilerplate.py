# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Copyright (c) 2020:
#     Laboratory for Computational Cellular Mechanobiology
#     Cuncheng Zhu (cuzhu@eng.ucsd.edu)
#     Christopher T. Lee (ctlee@ucsd.edu)
#     Ravi Ramamoorthi (ravir@cs.ucsd.edu)
#     Padmini Rangamani (prangamani@eng.ucsd.edu)
#

import pymem3dg as dg
import numpy as np
import numpy.typing as npt
import pymem3dg.util as dg_util


def getGeodesicDistance(
    face: npt.NDArray[np.int64], vertex: npt.NDArray[np.float64], point: list
):
    system = dg.System(face, vertex)
    system.parameters.point.pt = point
    if system.parameters.point.isFloatVertex:
        system.findFloatCenter()
    else:
        system.findVertexCenter()
    distance = system.computeGeodesicDistance()
    return distance


def prescribeGeodesicPoteinDensityDistribution(
    time: float,
    vertexMeanCuravtures: npt.NDArray[np.float64],
    geodesicDistance: npt.NDArray[np.float64],
):
    return dg_util.tanhDistribution(x=geodesicDistance, sharpness=10, center=1)
