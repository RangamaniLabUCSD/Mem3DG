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
#include <cstdarg>
#include <cstddef>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/src/Core/util/Constants.h"

#include "mem3dg/solver/geometry.h"
#include "mem3dg/solver/integrator/conjugate_gradient.h"
#include "mem3dg/solver/mesh_process.h"
#include "mem3dg/solver/system.h"
#include "mem3dg/version.h"

#include <geometrycentral/utilities/eigen_interop_helpers.h>

#include "mem3dg/mem3dg"
#include "pybind11/cast.h"

namespace gc = ::geometrycentral;
namespace mem3dg {
namespace solver {
namespace integrator {
namespace py = pybind11;

void init_parameters(py::module_ &pymem3dg) {
  // ==========================================================
  // =============   Simulation parameters      ===============
  // ==========================================================
  py::class_<Parameters::Boundary> boundary(pymem3dg, "Boundary", R"delim(
        The boundary conditions
    )delim");
  boundary.def(py::init<>());
  boundary.def_readwrite("shapeBoundaryCondition",
                         &Parameters::Boundary::shapeBoundaryCondition,
                         R"delim(
          get the option of "roller", "pin", "fixed", or "none" to specify shape boundary condition.
      )delim");
  boundary.def_readwrite("proteinBoundaryCondition",
                         &Parameters::Boundary::proteinBoundaryCondition,
                         R"delim(
          get the option of "pin", or "none" to specify protein boundary condition.
      )delim");

  py::class_<Parameters::Variation> variation(pymem3dg, "Variation", R"delim(
        Variation
    )delim");
  variation.def(py::init<>());
  variation.def_readwrite("isProteinVariation",
                          &Parameters::Variation::isProteinVariation,
                          R"delim(
          get the option of whether simulate protein variation
      )delim");
  variation.def_readwrite("isProteinConservation",
                          &Parameters::Variation::isProteinConservation,
                          R"delim(
          get the option of whether conserve protein mass
      )delim");
  variation.def_readwrite("isShapeVariation",
                          &Parameters::Variation::isShapeVariation,
                          R"delim(
          get the option of whether simulate shape variation
      )delim");
  variation.def_readwrite("geodesicMask", &Parameters::Variation::geodesicMask,
                          R"delim(
          get domain of shape variation
      )delim");
  variation.def_readwrite("updateMaskPeriod",
                          &Parameters::Variation::updateMaskPeriod,
                          R"delim(
          period of updating mask. measured in the unit of # of iterations
      )delim");

  py::class_<Parameters::Bending> bending(pymem3dg, "Bending", R"delim(
        The bending parameters
    )delim");
  bending.def(py::init<>());
  bending.def_readwrite("D", &Parameters::Bending::D,
                        R"delim(
          get thickness of the membrane (area difference model)
      )delim");
  bending.def_readwrite("alpha", &Parameters::Bending::alpha,
                        R"delim(
          get unity modulus of the membrane (area difference model)
      )delim");
  bending.def_readwrite("dA0", &Parameters::Bending::dA0,
                        R"delim(
          get preferred area difference of the membrane (area difference model)
      )delim");
  bending.def_readwrite("Kd", &Parameters::Bending::Kd,
                        R"delim(
          get deviatoric rigidity of the membrane
      )delim");
  bending.def_readwrite("Kdc", &Parameters::Bending::Kdc,
                        R"delim(
          get constant of deviatoric modulus vs protein density
      )delim");
  bending.def_readwrite("Kb", &Parameters::Bending::Kb,
                        R"delim(
          get Bending rigidity of the bare membrane
      )delim");
  bending.def_readwrite("Kbc", &Parameters::Bending::Kbc,
                        R"delim(
          get constant of bending modulus vs protein density
      )delim");
  bending.def_readwrite("H0c", &Parameters::Bending::H0c,
                        R"delim(
          get constant of spontaneous curvature vs protein density
      )delim");
  bending.def_readwrite("relation", &Parameters::Bending::relation,
                        R"delim(
          get relation between H0 and protein densit, "linear" or "hill"
      )delim");

  py::class_<Parameters::Tension> tension(pymem3dg, "Tension", R"delim(
        The surface tension parameters
    )delim");
  tension.def(py::init<>());
  tension.def_readwrite("form", &Parameters::Tension::form,
                        R"delim(
          functional to set the tension area relation
        args:
            total surface area of the mesh
        return:
            tuple of surface tension of the system and surface energy
      )delim");

  py::class_<Parameters::Osmotic> osmotic(pymem3dg, "Osmotic", R"delim(
        The osmotic pressure parameters
    )delim");
  osmotic.def_readwrite("form", &Parameters::Osmotic::form,
                        R"delim(
          functional to set the pressure volume relation
        args:
            enclosed volume of the mesh
        return:
            tuple of osmotic pressure of the system and pressure energy
      )delim");

  py::class_<Parameters::Adsorption> adsorption(pymem3dg, "Adsorption",
                                                R"delim(
        The adsorption parameters
    )delim");
  adsorption.def_readwrite("epsilon", &Parameters::Adsorption::epsilon,
                           R"delim(
          get adsorption energy per protein
      )delim");

  py::class_<Parameters::Aggregation> aggregation(pymem3dg, "Aggregation",
                                                  R"delim(
        The aggregation parameters
    )delim");
  aggregation.def_readwrite("chi", &Parameters::Aggregation::chi,
                            R"delim(
          get aggregation energy
      )delim");

  py::class_<Parameters::Entropy> entropy(pymem3dg, "Entropy",
                                          R"delim(
        The entropy parameters
    )delim");
  entropy.def_readwrite("xi", &Parameters::Entropy::xi,
                        R"delim(
          get entropy parameters
      )delim");

  py::class_<Parameters::External> external(pymem3dg, "External",
                                            R"delim(
        The external force parameters
    )delim");
  external.def_readwrite("form", &Parameters::External::form,
                         R"delim(
          functional to set the external force prescription
        args:
            vertexPositions
            vertexDualAreas
            time,
            geodesicDistance
        return:
            externalForce
      )delim");

  py::class_<Parameters::DPD> dpd(pymem3dg, "DPD",
                                  R"delim(
        The DPD parameters
    )delim");
  dpd.def_readwrite("gamma", &Parameters::DPD::gamma,
                    R"delim(
          get Dissipation coefficient
      )delim");

  py::class_<Parameters::Dirichlet> dirichlet(pymem3dg, "Dirichlet",
                                              R"delim(
        The Dirichlet energy parameters
    )delim");
  dirichlet.def_readwrite("eta", &Parameters::Dirichlet::eta,
                          R"delim(
          get coefficient
      )delim");

  py::class_<Parameters::SelfAvoidance> selfAvoidance(pymem3dg, "SelfAvoidance",
                                                      R"delim(
        The SelfAvoidance energy parameters
    )delim");
  selfAvoidance.def_readwrite("d", &Parameters::SelfAvoidance::d,
                              R"delim(
          get coefficient of limit distance
      )delim");
  selfAvoidance.def_readwrite("mu", &Parameters::SelfAvoidance::mu,
                              R"delim(
          get coefficient of penalty coefficient
      )delim");
  selfAvoidance.def_readwrite("n", &Parameters::SelfAvoidance::n,
                              R"delim(
          get the number excluding neighborhood layers
      )delim");
  selfAvoidance.def_readwrite("p", &Parameters::SelfAvoidance::p,
                              R"delim(
          get the period factor of self-avoidance computation
      )delim");

  py::class_<Parameters::Point> point(pymem3dg, "Point",
                                      R"delim(
        The Point energy parameters
    )delim");
  point.def_readwrite("prescribeNotableVertex",
                      &Parameters::Point::prescribeNotableVertex,
                      R"delim(
          functional to find the notable vertex of the mesh
        args:
            faceMatrix (npt.NDarray[int64])
            vertexMatrix (npt.NDarray[float64])
            geodesicDistance (list)
        return:
            notable vertex (bool list)
      )delim");
  point.def_readwrite("updateGeodesicsPeriod",
                      &Parameters::Point::updateGeodesicsPeriod, R"delim(
            the period factor of updating geodesics distance around notable vertices. Measured in the unit of # of iterations.
    )delim");
  point.def_readwrite("updateNotableVertexPeriod",
                      &Parameters::Point::updateNotableVertexPeriod, R"delim(
            the period factor of updating the vertex data of notable vertices. Measured in the unit of # of iterations.
    )delim");

  py::class_<Parameters::Protein> protein(pymem3dg, "Protein",
                                          R"delim(
        The protein distribution parameters
    )delim");
  protein.def_readwrite(
      "updateProteinDensityDistributionPeriod",
      &Parameters::Protein::updateProteinDensityDistributionPeriod,
      R"delim(
          period of updating protein density distribution. measured in # of iterations
      )delim");
  protein.def_readwrite("proteinInteriorPenalty",
                        &Parameters::Protein::proteinInteriorPenalty,
                        R"delim(
          get interior point parameter for protein density
      )delim");
  protein.def_readwrite(
      "prescribeProteinDensityDistribution",
      &Parameters::Protein::prescribeProteinDensityDistribution,
      R"delim(
          functional to set the protein density distribution prescription
        args:
            time (float)
            vertexMeanCurvatures (list)
            geodesicDistance (list)
        return:
            proteinDensity (list)
      )delim");

  py::class_<Parameters::Spring> spring(pymem3dg, "spring",
                                        R"delim(
        mesh spring forces parameters
    )delim");
  spring.def_readwrite("Kst", &Parameters::Spring::Kst,
                       R"delim(
          get Vertex shifting constant
      )delim");
  spring.def_readwrite("Ksl", &Parameters::Spring::Ksl,
                       R"delim(
          get Local stretching modulus
      )delim");
  spring.def_readwrite("Kse", &Parameters::Spring::Kse,
                       R"delim(
          get Edge spring constant
      )delim");

  py::class_<Parameters> parameters(pymem3dg, "Parameters", R"delim(
        The parameters
    )delim");
  parameters.def(py::init<>());
  parameters.def_readwrite("bending", &Parameters::bending,
                           R"delim(
          bending parameters
      )delim");
  parameters.def_readwrite("tension", &Parameters::tension,
                           R"delim(
          tension parameters
      )delim");
  parameters.def_readwrite("osmotic", &Parameters::osmotic,
                           R"delim(
        "osmotic parameters
      )delim");
  parameters.def_readwrite("adsorption", &Parameters::adsorption,
                           R"delim(
          adsorption parameters
      )delim");
  parameters.def_readwrite("aggregation", &Parameters::aggregation,
                           R"delim(
          aggregation parameters
      )delim");
  parameters.def_readwrite("entropy", &Parameters::entropy,
                           R"delim(
          entropy parameters
      )delim");
  parameters.def_readwrite("dirichlet", &Parameters::dirichlet,
                           R"delim(
          dirichlet parameters
      )delim");
  parameters.def_readwrite("selfAvoidance", &Parameters::selfAvoidance,
                           R"delim(
          selfAvoidance parameters
      )delim");
  parameters.def_readwrite("dpd", &Parameters::dpd,
                           R"delim(
          dpd parameters
      )delim");
  parameters.def_readwrite("external", &Parameters::external,
                           R"delim(
          external parameters
      )delim");
  parameters.def_readwrite("boundary", &Parameters::boundary,
                           R"delim(
          boundary parameters
      )delim");
  parameters.def_readwrite("point", &Parameters::point,
                           R"delim(
          point parameters
      )delim");
  parameters.def_readwrite("protein", &Parameters::protein,
                           R"delim(
          protein parameters
      )delim");
  parameters.def_readwrite("variation", &Parameters::variation,
                           R"delim(
          variation parameters
      )delim");
  parameters.def_readwrite("temperature", &Parameters::temperature,
                           R"delim(
          get Temperature
      )delim");
  parameters.def_readwrite("proteinMobility", &Parameters::proteinMobility,
                           R"delim(
          get protein mobility constant
      )delim");
  parameters.def_readwrite("damping", &Parameters::damping,
                           R"delim(
          get damping constant
      )delim");
  parameters.def_readwrite("spring", &Parameters::spring,
                           R"delim(
          get spring parameters
      )delim");
}
} // namespace integrator
} // namespace solver
} // namespace mem3dg
