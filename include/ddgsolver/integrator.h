#pragma once

#include <geometrycentral/surface/geometry.h>
#include <geometrycentral/surface/rich_surface_mesh_data.h>

#include "force.h"

namespace ddgsolver {
  namespace integration {

    DLL_PUBLIC double getTotalEnergy(Force& f);

    DLL_PUBLIC void stormerVerlet(Force& f, double dt, double total_time,
      double tolerance);

    DLL_PUBLIC void velocityVerlet(Force& f, double dt, double total_time,
      double tolerance, double tSave, std::string outputDir);

    DLL_PUBLIC void getLogFiles(Force& f, double dt, double total_time,
      double tolerance, double tSave, std::string outputDir);

  }// namespace integration
} // namespace ddgsolver
