#pragma once

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "force.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;


class DLL_PUBLIC integrator {
public:
  double &timeStep;
  double &timeSpan;
  gcs::HalfedgeMesh &mesh;
  gcs::VertexPositionGeometry &vpg;
  Parameters &p;
  ddgsolver::Force &f;
  double tolerance;
  double bendingEnergy;
  double pastBendingEnergy;

  integrator(gcs::HalfedgeMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
             ddgsolver::Force &f_, double &h, double &T, Parameters &p_,
             double eps)
      : mesh(mesh_), vpg(vpg_), f(f_), timeStep(h), timeSpan(T), p(p_),
        tolerance(eps) {
    bendingEnergy = 0;
    pastBendingEnergy = 0;
  }

  void stormerVerlet();
  void velocityVerlet();
  void getBendingEnergy();
};

DLL_PUBLIC void velocityVerlet(Force &f, double dt, double total_time,
                               double tolerance);

} // namespace ddgsolver
