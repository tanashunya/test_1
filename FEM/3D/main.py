import MOTOR_FEM
import os

if not os.path.exists("./vtk"):
    os.mkdir("./vtk")

obj = MOTOR_FEM.MOTOR_FEM()
obj.main_func()