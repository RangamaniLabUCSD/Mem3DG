#pragma once

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "force.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

struct DLL_PUBLIC parameters {
  double Kb;    /// Bending modulus
  double H0;    /// Spontaneous curvature
  double Ksl;   /// Local stretching modulus
  double Ksg;   /// Global stretching modulus
  double Kse;   /// 
  double Kv;    /// Volume regularization
  double gamma; /// Dissipation coefficient
  double Vt;    /// Reduced volume
  double kt;    /// Boltzmann constant*Temperature
  double sigma; /// Noise
};

class DLL_PUBLIC integrator {
public:
  double &timeStep;
  double &timeSpan;
  gcs::HalfedgeMesh &mesh;
  gcs::VertexPositionGeometry &vpg;
  parameters &p;
  ddgsolver::Force &f;
  double tolerance;
  double bendingEnergy;
  double pastBendingEnergy;

  integrator(gcs::HalfedgeMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
             ddgsolver::Force &f_, double &h, double &T, parameters &p_,
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
} // namespace ddgsolver
