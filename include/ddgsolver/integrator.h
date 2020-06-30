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
  double &dt;
  double &total_time;
  gcs::HalfedgeMesh &mesh;
  gcs::VertexPositionGeometry &vpg;
  gcs::PlyHalfedgeMeshData &plyData;
  Parameters &p;
  ddgsolver::Force &f;
  double tolerance;
  double totalEnergy;
  double pastTotalEnergy;
  double tSave;

  integrator(gcs::HalfedgeMesh &mesh_, gcs::VertexPositionGeometry &vpg_,
    gcs::PlyHalfedgeMeshData& plyData_, ddgsolver::Force &f_, double &h, 
    double &T, Parameters &p_, double eps, double tSave_)
      : mesh(mesh_), vpg(vpg_), plyData(plyData_), f(f_), dt(h), total_time(T),
    p(p_), tolerance(eps), tSave(tSave_) {
    totalEnergy = 0;
    pastTotalEnergy = 0;
  }

  void stormerVerlet();
  void velocityVerlet();
  void getTotalEnergy();
  void getLogFiles();

};

DLL_PUBLIC void velocityVerlet(Force &f, double dt, double total_time,
                               double tolerance, double tSave = 0.2);

} // namespace ddgsolver
