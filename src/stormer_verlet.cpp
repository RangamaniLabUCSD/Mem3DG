#include "ddgsolver/integrator.h"
#include "ddgsolver/force.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/vector3.h>

#include <iostream>

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

	void integrator::stormerVerlet() {
		gcs::FaceData<size_t> faceInd = vpg.faceIndices;
		gc::Vector3 totalForce;
		for (size_t i = 0; i < timeSpan / timeStep; i++) {
			/*polyscope::registerSurfaceMesh("myMesh",
				ptrvpg->inputVertexPositions,
				ptrmesh->getFaceVertexList());*/
				//polyscope::show();
			f.getVelocityFromPastPosition(timeStep);
			f.getBendingForces();
			f.getStretchingForces();
			f.getPressureForces();
			f.getDampingForces();
			f.getStochasticForces();

			gcs::VertexData<gc::Vector3> temp = vpg.inputVertexPositions;
			for (gcs::Vertex v : mesh.vertices()) {
				bool flag = true;
				for (gcs::Face f : v.adjacentFaces()) {
					if (faceInd[f] == 0) {
						flag = true; // change it to false to have time integration excluded for a facet of vertices 
					}
				}
				if (flag == true) {
					vpg.inputVertexPositions[v] *= 2;
					totalForce = f.bendingForces[v]
						+ f.stretchingForces[v]
						+ f.pressureForces[v]
						+ f.dampingForces[v]
						+ f.stochasticForces[v];
					vpg.inputVertexPositions[v] += totalForce * timeStep * timeStep - f.pastPositions[v];
				}
			}
			//std::cout << "total force:  " << totalForce.norm() << std::endl;
			f.update_Vertex_positions();
			f.pastPositions = temp;
			getBendingEnergy();
			if (abs(pastBendingEnergy - bendingEnergy) / bendingEnergy < 1e-6) { break; }
			std::cout << "energy: " << bendingEnergy << std::endl;
			std::cout << "process: " << int(double(i) / (timeSpan / timeStep) * 100) << "%" << std::endl;
		}
	}
}
