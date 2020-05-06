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

	void integrator::velocityVerlet() {
		f.vertexVelocity.fill({ 0,0,0 });
		gcs::VertexData<gc::Vector3> totalForce(mesh,{ 0,0,0 });
		gcs::VertexData<gc::Vector3> newTotalForce(mesh,{ 0,0,0 });
		gc::Vector3 COMVelocity;
		Eigen::Map<Eigen::Matrix<double, 1, 3>> COMVelocity_e(&COMVelocity[0]);
		//auto COMVelocity_e = ddgsolver::EigenMap<double, 3>(COMVelocity);
		for (size_t i = 0; i < timeSpan / timeStep; i++) {

			f.getBendingForces(p.Kb, p.H0);
			f.getStretchingForces(p.Ksl, p.Ksg);
			f.getPressureForces(p.Kv, p.Vt);
			f.getDampingForces(p.gamma);
			f.getStochasticForces(p.sigma);

			COMVelocity_e = ddgsolver::EigenMap<double, 3>(f.vertexVelocity).colwise().sum() / mesh.nVertices();
			for (gcs::Vertex v : mesh.vertices()) {
				vpg.inputVertexPositions[v] += (f.vertexVelocity[v] - COMVelocity) * timeStep
					+ totalForce[v] * timeStep * timeStep * 0.5;

				newTotalForce[v] = f.bendingForces[v]
					+ f.stretchingForces[v]
					+ f.pressureForces[v]
					+ f.dampingForces[v]
					+ f.stochasticForces[v];

				f.vertexVelocity[v] += (totalForce[v] + newTotalForce[v]) * timeStep * 0.5;
				totalForce[v] = newTotalForce[v];
			}
			f.update_Vertex_positions();
			/*double energy = getBendingEnergy();
			std::cout << energy << std::endl;*/
		}
	}
}