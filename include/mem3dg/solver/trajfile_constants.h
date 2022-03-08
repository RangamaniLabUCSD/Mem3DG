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

/**
 * @file  mutable_trajfile.h
 * @brief Netcdf trajectory output support
 *
 */

#pragma once

namespace mem3dg {
namespace solver {

#ifdef MEM3DG_WITH_NETCDF
// DIMENSIONS NAMES
static const std::string POLYGON_ORDER_NAME = "polygon_dims";
/// Number of vertices per polygon
static const std::size_t POLYGON_ORDER = 3;

static const std::string NPOLYGONS_NAME = "npolygons";

static const std::string SPATIAL_DIMS_NAME = "spatial";
/// Number of spatial dimensions
static const std::size_t SPATIAL_DIMS = 3;
/// Name of frames
static const std::string FRAME_NAME = "frame";
/// nvertices
static const std::string NVERTICES_NAME = "nvertices";
/// nconers
static const std::string NCORNERS_NAME = "ncorners";

/// Name of conventions
static const std::string CONVENTIONS_NAME = "Conventions";
/// Convention value
static const std::string CONVENTIONS_VALUE = "Mem3DG";
/// Conventions version
static const std::string CONVENTIONS_VERSION_NAME = "ConventionsVersion";
/// Conventions version value
static const std::string CONVENTIONS_VERSION_VALUE = "0.0.1";

static const std::string PARAM_GROUP_NAME = "Parameters";
static const std::string TRAJ_GROUP_NAME = "Trajectory";

/// Name of the units labels
static const std::string UNITS = "units";
/// Value for length units
static const std::string LEN_UNITS = " micrometers ";
/// Value for time units
static const std::string TIME_UNITS = " seconds ";
/// Value for force units
static const std::string FORCE_UNITS = " nanonewtons ";

// Data/Variable block names
/// Name of time data
static const std::string TIME_VAR = "time";
/// Name of isSmooth data
static const std::string ISSMOOTH_VAR = "issmooth";
/// Name of coordinates data
static const std::string COORD_VAR = "coordinates";
/// Name of the mesh topology data
static const std::string TOPO_VAR = "topology";
/// Name of the mesh topology data
static const std::string TOPO_FRAME_VAR = "topologyframe";
/// Name of the mesh corner angle data
static const std::string ANGLE_VAR = "angle";
/// Name of the refMesh coordinates data
static const std::string REFCOORD_VAR = "refcoordinates";
/// Name of the velocity data
static const std::string VEL_VAR = "velocities";
/// Name of the velocity data
static const std::string EXTF_VAR = "externalForce";
/// Name of the mean curvature data
static const std::string MEANCURVE_VAR = "meancurvature";
/// Name of the Gaussian curvature data
static const std::string GAUSSCURVE_VAR = "gausscurvature";
/// Name of the protein density data
static const std::string PHI_VAR = "proteindensity";
/// Name of the spontaneous curvature data
static const std::string SPONCURVE_VAR = "sponcurvature";
/// Name of the external pressure data
static const std::string EXTERNFORCE_VAR = "externpressure";
/// Name of the chemical potential data
static const std::string CHEMPOTENTIAL_VAR = "chempotential";
/// Name of the physical pressure data
static const std::string PHYSFORCE_VAR = "physpressure";
/// Name of the capillary pressure data
static const std::string CAPFORCE_VAR = "cappressure";
/// Name of the bending pressure data
static const std::string BENDFORCE_VAR = "bendpressure";
/// Name of the inside pressure data
static const std::string OSMOTICFORCE_VAR = "insidepressure";
/// Name of the line tension pressure data
static const std::string LINEFORCE_VAR = "linepressure";
/// Name of the bending energy data
static const std::string BENDENER_VAR = "bendenergy";
/// Name of the surface energy data
static const std::string SURFENER_VAR = "surfenergy";
/// Name of the pressure energy data
static const std::string PRESSENER_VAR = "pressenergy";
/// Name of the kinetic energy data
static const std::string KINEENER_VAR = "kineenergy";
/// Name of the adsorption energy data
static const std::string ADSPENER_VAR = "adspenergy";
/// Name of the line tension energy data
static const std::string LINEENER_VAR = "lineenergy";
/// Name of the chemical energy data
static const std::string TOTALENER_VAR = "totalenergy";
/// Name of the  chem Error Norm data
static const std::string CHEMERRORNORM_VAR = "chemerrornorm";
/// Name of the  Error Norm data
static const std::string ERRORNORM_VAR = "errornorm";
/// Name of the  Error Norm data
static const std::string BENDNORM_VAR = "bendnorm";
/// Name of the  Error Norm data
static const std::string SURFNORM_VAR = "surfnorm";
/// Name of the  Error Norm data
static const std::string PRESSNORM_VAR = "pressnorm";
/// Name of the  Error Norm data
static const std::string LINENORM_VAR = "linenorm";
/// Name of the volume data
static const std::string VOLUME_VAR = "volume";
/// Name of the height data
static const std::string HEIGHT_VAR = "height";
/// Name of the surface area data
static const std::string SURFAREA_VAR = "surfacearea";
/// Name of the reference volume data
static const std::string REFVOLUME_VAR = "refvolume";
/// Name of the reference surface area data
static const std::string REFSURFAREA_VAR = "refsurfarea";
/// Name of the mask data
static const std::string MASK_VAR = "mask";
/// Name of the curvature difference data
static const std::string H_H0_VAR = "curvaturediff";

/// Name of uint array vlen type
static const std::string UINT_ARR = "uint_array";
/// Name of double array vlen type
static const std::string DOUBLE_ARR = "double_array";

#endif
} // namespace solver
} // namespace mem3dg
