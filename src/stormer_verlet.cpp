#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <iostream>

namespace ddgsolver {
	namespace integration {
		namespace gc = ::geometrycentral;
		namespace gcs = ::geometrycentral::surface;

		void stormerVerlet(Force& f, double dt, double total_time,
			double tolerance) {
			gcs::FaceData<size_t> faceInd = f.vpg.faceIndices;
			gc::Vector3 totalForce;
			for (size_t i = 0; i < total_time / dt; i++) {
				/*polyscope::registerSurfaceMesh("myMesh",
					ptrVpg->inputVertexPositions,
					ptrMesh->getFaceVertexList());*/
					//polyscope::show();
				f.getVelocityFromPastPosition(dt);
				f.getBendingForces();
				f.getStretchingForces();
				f.getPressureForces();
				f.getDPDForces();
				f.getExternalForces();

				gcs::VertexData<gc::Vector3> temp = f.vpg.inputVertexPositions;
				for (gcs::Vertex v : f.mesh.vertices()) {
					bool flag = true;
					for (gcs::Face f : v.adjacentFaces()) {
						if (faceInd[f] == 0) {
							flag = true; // change it to false to have time integration excluded for a facet of vertices 
						}
					}
					if (flag == true) {
						f.vpg.inputVertexPositions[v] *= 2;
						totalForce = f.bendingForces[v]
							+ f.stretchingForces[v]
							+ f.pressureForces[v]
							+ f.dampingForces[v]
							+ f.stochasticForces[v]
							+ f.externalForces[v];
						f.vpg.inputVertexPositions[v] += totalForce * dt * dt - f.pastPositions[v];
					}
				}
				//std::cout << "total force:  " << totalForce.norm() << std::endl;
				f.update_Vertex_positions();
				f.pastPositions = temp;
				double totalEnergy = getTotalEnergy(f);
				
				std::cout << "energy: " << totalEnergy << std::endl;
				std::cout << "process: " << int(double(i) / (total_time / dt) * 100) << "%" << std::endl;
			}
		}
	} // namespace integration
} // namespace ddgsolver
