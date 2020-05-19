#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <Eigen/Core>

#include <pcg_random.hpp>

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

		std::default_random_engine random_generator;
		std::uniform_int_distribution<int> uniform_dist(1, 10);

		for (size_t i = 0; i < timeSpan / timeStep; i++) {

			double timeStepHere = timeStep * (1 / pow(i, 1 / 2));
			p.sigma = sqrt(2 * p.gamma * p.kt / timeStepHere);

			int prob = 4;

			int flag = uniform_dist(random_generator);
			if (true){ f.getBendingForces(p.Kb, p.H0); }

			flag = uniform_dist(random_generator);
			if (flag < prob) { f.getStretchingForces(p.Ksl, p.Ksg, p.Kse); }
			else { f.stretchingForces.fill({ 0.0,0.0,0.0 }); }
			
			flag = uniform_dist(random_generator);
			if (flag < prob) { f.getPressureForces(p.Kv, p.Vt); }
			else { f.pressureForces.fill({ 0.0,0.0,0.0 }); }
			
			flag = uniform_dist(random_generator);
			if (flag < prob) {
				f.getDampingForces(p.gamma);
				f.getStochasticForces(p.sigma);
			}
			else { 
				f.dampingForces.fill({ 0.0,0.0,0.0 }); 
				f.stochasticForces.fill({ 0.0,0.0,0.0 });
			}
			//f.getBendingForces(p.Kb, p.H0);
			//f.getStretchingForces(p.Ksl, p.Ksg,p.Kse);
			//f.getPressureForces(p.Kv, p.Vt);
			//f.getDampingForces(p.gamma);
			//f.getStochasticForces(p.sigma);

			COMVelocity_e = ddgsolver::EigenMap<double, 3>(f.vertexVelocity).colwise().sum() / mesh.nVertices();
			for (gcs::Vertex v : mesh.vertices()) {
				vpg.inputVertexPositions[v] += (f.vertexVelocity[v] - COMVelocity) * timeStepHere
					+ totalForce[v] * timeStepHere * timeStepHere * 0.5;

				newTotalForce[v] = f.bendingForces[v]
					+ f.stretchingForces[v]
					+ f.pressureForces[v]
					+ f.dampingForces[v]
					+ f.stochasticForces[v];

				/*std::cout << "bf: " << f.bendingForces[v].norm()
					<< "sf: " << f.stretchingForces[v].norm()
					<< "pf: " << f.pressureForces[v].norm()
					<< "df: " << f.dampingForces[v].norm()
					<< "xf: " << f.stochasticForces[v].norm() <<std::endl;*/

				f.vertexVelocity[v] += (totalForce[v] + newTotalForce[v]) * timeStepHere * 0.5;
				totalForce[v] = newTotalForce[v];
			}
			f.update_Vertex_positions();

			getBendingEnergy();
			if (((abs(pastBendingEnergy - bendingEnergy) / bendingEnergy) < tolerance)
				&& (i>1) && (abs(f.volume - f.targetVolume * p.Vt) / (f.targetVolume * p.Vt) < 1e-2)
				&& (abs(f.surfaceArea - f.targetSurfaceArea) / (f.targetSurfaceArea) < 1e-2)) { break; }
			// std::cout << "energy: " << bendingEnergy << std::endl;

			//std::cout << "process: " << int(double(i) / (timeSpan / timeStep) * 100) << "%" << std::endl;
			std::cout << "process: " << i << std::endl;

		}
	}
}