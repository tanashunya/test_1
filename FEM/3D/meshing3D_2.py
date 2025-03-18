import pygmsh
import gmsh
from meshing3D import visualize_mesh_by_pyvista


# Clear previous model
mesh_size = 0.1
geom = pygmsh.occ.Geometry()
model3D = geom.__enter__()
box0 = model3D.add_box([0.0, 0, 0], [1, 1, 1])
box1 = model3D.add_box([0.5, 0.5, 1], [0.5, 0.5, 1])
ball = model3D.add_ball([0.5, 0.5, 0.5], 0.25)

union = model3D.boolean_union([box0, box1])
union_minus_ball = model3D.boolean_fragments(union, ball)
model3D.synchronize()

model3D.add_physical(union, "Union")
model3D.add_physical(union_minus_ball, "Union minus ball")

geom.generate_mesh(dim=3)
gmsh.write("mesh3D.msh")
model3D.__exit__()

visualize_mesh_by_pyvista("mesh3D.msh")