# 3次元の立方体領域を構造格子で切る
import gmsh

gmsh.initialize()
gmsh.model.add("cube")

# Points
points = [
    (0, 0, 0, 1.0),
    (1, 0, 0, 1.0),
    (1, 1, 0, 1.0),
    (0, 1, 0, 1.0),
    (0, 0, 1, 1.0),
    (1, 0, 1, 1.0),
    (1, 1, 1, 1.0),
    (0, 1, 1, 1.0)
]
point_tags = [gmsh.model.geo.addPoint(*p) for p in points]

# Lines
lines = [
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 1),
    (1, 5),
    (2, 6),
    (3, 7),
    (4, 8),
    (8, 7),
    (7, 6),
    (6, 5),
    (5, 8)
]
line_tags = [gmsh.model.geo.addLine(*l) for l in lines]

# Line Loops and Surfaces
line_loops = [
    [1, 2, 3, 4],
    [3, 8, 9, -7],
    [4, 5, 12, -8],
    [7, 10, -6, 2],
    [1, 6, 11, -5],
    [9, 10, 11, 12]
]
surface_tags = [gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(ll)]) for ll in line_loops]

# Surface Loop and Volume
surface_loop_tag = gmsh.model.geo.addSurfaceLoop(surface_tags)
volume_tag = gmsh.model.geo.addVolume([surface_loop_tag])

# Transfinite Lines
for tag in line_tags:
    gmsh.model.geo.mesh.setTransfiniteCurve(tag, 21, coef=1.0)

# Transfinite Surfaces
for tag in surface_tags:
    gmsh.model.geo.mesh.setTransfiniteSurface(tag)

# Transfinite Volume
gmsh.model.geo.mesh.setTransfiniteVolume(volume_tag)

# Physical Groups
gmsh.model.addPhysicalGroup(2, [surface_tags[5]], name="movingWall")
gmsh.model.addPhysicalGroup(2, surface_tags[:5], name="wall")
gmsh.model.addPhysicalGroup(3, [volume_tag], name="fluid")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)

# Save and finalize
gmsh.write("cube4.msh1")
gmsh.finalize()