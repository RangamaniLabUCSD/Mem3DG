#include "ddgsolver/integrator.h"

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include "force.h"

namespace ddgsolver {
namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

void stormerVerlet(ddgsolver::Force& f, ddgsolver::parameters& p) {
	for (size_t i = 0; i < T / h; i++) {
		/*polyscope::registerSurfaceMesh("myMesh",
			ptrvpg->inputVertexPositions,
			ptrmesh->getFaceVertexList());*/
			//polyscope::show();
		f.getBendingForces(Kb, H0);
		f.getStretchingForces(Ksl, Ksg);
		f.getPressureForces(Kv, Vt);
		f.getDampingForces(gamma);
		f.getStochasticForces(sigma);

		gcs::VertexData<gc::Vector3> temp = vpg.inputVertexPositions;
		for (gcs::Vertex v : mesh.vertices()) {
			vpg.inputVertexPositions[v] *= 2;
			vpg.inputVertexPositions[v] += (f.bendingForces[v]
				+ f.stretchingForces[v]
				+ f.pressureForces[v]
				+ f.dampingForces[v]
				+ f.stochasticForces[v]) * h * h - f.pastPositions[v];
		}
		f.update_Vertex_positions();
		f.pastPositions = temp;
}
}
