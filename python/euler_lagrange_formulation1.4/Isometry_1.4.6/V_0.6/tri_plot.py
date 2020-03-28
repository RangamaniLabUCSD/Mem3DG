import trimesh
import numpy as np

import threevis
import openmesh as om



v = np.load("vertices.npy")
f = np.load("faces.npy")
# mesh = trimesh.Trimesh(vertices=v,
#                         faces=f)
#mesh.show()



# m = om.read_trimesh("v_06.obj")
# tv.display_openmesh(m)

#protverts, protedges, protfaces = meshes[0].to_ndarray()
#bverts, bedges, bfaces = meshes[1].to_ndarray()
ctx = threevis.Context(width=640, height=480)
colors = np.zeros((len(f), 3))
colors[:,0] = np.ones(len(f))
ctx.draw_faces(v,f, colors=threevis.FaceAttribute(colors))
#ctx.draw_edges(protverts, protedges)
ctx.display()