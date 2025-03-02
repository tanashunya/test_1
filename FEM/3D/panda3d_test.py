import vtk
from panda3d.core import GeomVertexFormat, GeomVertexData, Geom, GeomNode, GeomTriangles, GeomVertexWriter, NodePath
from direct.showbase.ShowBase import ShowBase

class VTKtoOBJConverter:
    def __init__(self, vtk_filename, obj_filename):
        self.vtk_filename = vtk_filename
        self.obj_filename = obj_filename

    def convert(self):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.vtk_filename)
        reader.Update()

        polydata = reader.GetOutput()

        with open(self.obj_filename, 'w') as obj_file:
            for i in range(polydata.GetNumberOfPoints()):
                point = polydata.GetPoint(i)
                obj_file.write(f'v {point[0]} {point[1]} {point[2]}\n')

            for i in range(polydata.GetNumberOfCells()):
                cell = polydata.GetCell(i)
                if cell.GetNumberOfPoints() == 3:
                    ids = cell.GetPointIds()
                    obj_file.write(f'f {ids.GetId(0) + 1} {ids.GetId(1) + 1} {ids.GetId(2) + 1}\n')

class Panda3DApp(ShowBase):
    def __init__(self, obj_filename):
        ShowBase.__init__(self)
        self.obj_filename = obj_filename
        self.load_obj()

    def load_obj(self):
        format = GeomVertexFormat.getV3()
        vdata = GeomVertexData('name', format, Geom.UHStatic)
        vertex = GeomVertexWriter(vdata, 'vertex')

        geom = Geom(vdata)
        tris = GeomTriangles(Geom.UHStatic)

        with open(self.obj_filename, 'r') as obj_file:
            vertices = []
            for line in obj_file:
                if line.startswith('v '):
                    parts = line.split()
                    vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
                elif line.startswith('f '):
                    parts = line.split()
                    tris.addVertices(int(parts[1]) - 1, int(parts[2]) - 1, int(parts[3]) - 1)

            for vertex_data in vertices:
                vertex.addData3f(*vertex_data)

        geom.addPrimitive(tris)
        node = GeomNode('gnode')
        node.addGeom(geom)
        node_path = NodePath(node)
        node_path.reparentTo(self.render)

if __name__ == '__main__':
    vtk_filename = 'vtk/material.vtk'
    obj_filename = 'magnet.obj'

    import vtk

    # # VTKファイルを作成する関数
    # def create_vtk_file(filename):
    #     # ポリゴンデータを作成
    #     points = vtk.vtkPoints()
    #     polygons = vtk.vtkCellArray()

    #     # 点を追加
    #     points.InsertNextPoint(0.0, 0.0, 0.0)
    #     points.InsertNextPoint(1.0, 0.0, 0.0)
    #     points.InsertNextPoint(1.0, 1.0, 0.0)
    #     points.InsertNextPoint(0.0, 1.0, 0.0)

    #     # ポリゴンを追加
    #     polygon = vtk.vtkPolygon()
    #     polygon.GetPointIds().SetNumberOfIds(4)
    #     polygon.GetPointIds().SetId(0, 0)
    #     polygon.GetPointIds().SetId(1, 1)
    #     polygon.GetPointIds().SetId(2, 2)
    #     polygon.GetPointIds().SetId(3, 3)
    #     polygons.InsertNextCell(polygon)

    #     # ポリゴンデータをポリデータに設定
    #     polydata = vtk.vtkPolyData()
    #     polydata.SetPoints(points)
    #     polydata.SetPolys(polygons)

    #     # VTKファイルに書き込む
    #     writer = vtk.vtkPolyDataWriter()
    #     writer.SetFileName(filename)
    #     writer.SetInputData(polydata)
    #     writer.Write()

    # # VTKファイルを作成
    # create_vtk_file(vtk_filename)

    # import pyvista as pv

    # # VTKファイルを読み込む
    # mesh = pv.read(vtk_filename)

    # # プロットする
    # plotter = pv.Plotter()
    # plotter.add_mesh(mesh, color='white')
    # plotter.show()

    converter = VTKtoOBJConverter(vtk_filename, obj_filename)
    converter.convert()

    app = Panda3DApp(obj_filename)
    app.run()
