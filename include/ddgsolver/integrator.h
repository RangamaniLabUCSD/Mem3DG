#pragma once

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include "force.h"

namespace ddgsolver {
	namespace gc = ::geometrycentral;
	namespace gcs = ::geometrycentral::surface;

	struct parameters {
		double Kb, H0, Ksl, Ksg, Kv, gamma, Vt, kt, sigma;
	};

	class integrator {
	public:
		double& timeStep;
		double& timeSpan;
		gcs::HalfedgeMesh& mesh;
		gcs::VertexPositionGeometry& vpg;
		parameters p;

		integrator(gcs::HalfedgeMesh& mesh_, gcs::VertexPositionGeometry& vpg_,
			double& h, double& T, parameters& p_) :mesh(mesh_), vpg(vpg_), timeStep(h), timeSpan(T), p(p_) {}

		void stormerVerlet(ddgsolver::Force& f);
		void velocityVerlet(ddgsolver::Force& f);

	};
}

