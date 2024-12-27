import FEM_3D
import os

if not os.path.exists("./vtk"):
    os.mkdir("./vtk")

obj = FEM_3D.FEM_3D()
obj.main_func()
