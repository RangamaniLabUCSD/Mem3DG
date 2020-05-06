#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <iostream>


namespace ddgsolver {
	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;
	double integrator::getBendingEnergy(double H0) { double Eb = 0; return Eb; }

	double integrator::getBendingEnergy() {
		double Eb = 0;
		auto positions = ddgsolver::EigenMap<double, 3>(vpg.inputVertexPositions);
		auto A = f.L.transpose() * f.M_inv * f.L; //could further cached since the same for all time pt
		for (size_t i = 0; i < positions.rows(); i++) {
			for (size_t j = 0; j < positions.rows(); j++) {
				Eb += positions.row(i) * A * positions.transpose().col(j);
			}
		}
		return Eb;
	}
}

