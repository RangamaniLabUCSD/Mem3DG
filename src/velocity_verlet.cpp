#include "ddgsolver/force.h"
#include "ddgsolver/integrator.h"
#include "ddgsolver/meshops.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>
#include <geometrycentral/surface/meshio.h>

#include <Eigen/Core>

#include <pcg_random.hpp>

#include <iostream>

namespace ddgsolver {
  namespace integration {
    namespace gc = ::geometrycentral;
    namespace gcs = ::geometrycentral::surface;

    void velocityVerlet(Force& f, double dt, double total_time,
      double tolerance, double tSave, std::string outputDir) {

      getLogFiles(f, dt, total_time, tolerance, tSave, outputDir);

      Eigen::Matrix<double, Eigen::Dynamic, 3> force;
      Eigen::Matrix<double, Eigen::Dynamic, 3> newForce;
      force.resize(f.mesh.nVertices(), 3);
      force.setZero();
      newForce.resize(f.mesh.nVertices(), 3);
      newForce.setZero();

      int nSave = int(tSave / dt);

      auto vel_e = EigenMap<double, 3>(f.vel);
      auto pos_e = EigenMap<double, 3>(f.vpg.inputVertexPositions);

      const double hdt = 0.5 * dt;
      const double hdt2 = hdt * dt;

      Eigen::Matrix<double, Eigen::Dynamic, 3> staticForce;
      Eigen::Matrix<double, Eigen::Dynamic, 3> dynamicForce;

      for (int i = 0; i <= total_time / dt; i++) {
        // Update all forces
        //f.getBendingForces();
        //f.getStretchingForces();
        //f.getPressureForces();
        f.getConservativeForces();
        f.getDPDForces();
        f.getExternalForces();

        //std::cout << "bf: " << EigenMap<double, 3>(f.bendingForces).norm()
        //  << "sf: " << EigenMap<double, 3>(f.stretchingForces).norm()
        //  << "pf: " << EigenMap<double, 3>(f.pressureForces).norm()
        //  << "df: " << EigenMap<double, 3>(f.dampingForces).norm()
        //  << "xf: " << EigenMap<double, 3>(f.stochasticForces).norm() << std::endl;

        pos_e +=
          (vel_e.rowwise() - (vel_e.colwise().sum() / f.mesh.nVertices())) * dt +
          force * hdt2;

        staticForce = EigenMap<double, 3>(f.bendingForces) +
          EigenMap<double, 3>(f.stretchingForces) +
          EigenMap<double, 3>(f.pressureForces) +
          EigenMap<double, 3>(f.externalForces);
        dynamicForce = EigenMap<double, 3>(f.dampingForces) +
          EigenMap<double, 3>(f.stochasticForces);
        newForce = staticForce + dynamicForce;

        vel_e += (force + newForce) * hdt;
        force = newForce;
        f.update_Vertex_positions(); // recompute cached values;
        double staticForce_mag = staticForce.norm();

        if (staticForce_mag < tolerance) {
          break;
        }

        if ((i % nSave == 0) || (i == int(total_time / dt))) {

          f.richData.addGeometry(f.vpg);

          gcs::VertexData<double> H(f.mesh);
          H.fromVector(f.M_inv * f.H);
          f.richData.addVertexProperty("mean_curvature", H);

          gcs::VertexData<double> f_ext(f.mesh);
          f_ext.fromVector(f.appliedForceMagnitude);
          f.richData.addVertexProperty("external_force", f_ext);

          gcs::VertexData<double> fn(f.mesh);
          fn.fromVector(rowwiseDotProduct(staticForce,
            EigenMap<double, 3>(f.vpg.vertexNormals)));
          f.richData.addVertexProperty("normal_force", fn);

          gcs::VertexData<double> ft(f.mesh);
          ft.fromVector((staticForce - rowwiseScaling(rowwiseDotProduct(staticForce,
            EigenMap<double, 3>(f.vpg.vertexNormals)),
            EigenMap<double, 3>(f.vpg.vertexNormals))).rowwise().norm());
          f.richData.addVertexProperty("tangential_force", ft);

          char buffer[50];
          sprintf(buffer, "t=%d.ply", int(i * dt * 100));
          f.richData.write(outputDir + buffer);

          std::cout << "time: " << i * dt << std::endl;
          std::cout << "force: " << staticForce_mag << std::endl;
          std::cout << "area: " << f.surfaceArea / f.initialSurfaceArea << std::endl;
          std::cout << "total volume:  " << f.volume / f.maxVolume / f.P.Vt << std::endl;

        }

      }

    }
  }// namespace integration
} // namespace ddgsolver
