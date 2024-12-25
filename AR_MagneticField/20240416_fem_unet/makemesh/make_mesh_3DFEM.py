import gmsh
import meshio

def create_3d_mesh(output_filename='mesh.msh'):
    gmsh.initialize()
    gmsh.model.add("3D Mesh")

    # ここにジオメトリの定義を追加します。
    # 例: 立方体の作成
    lc = 1.0  # メッシュの大きさを定義
    gmsh.model.occ.addBox(0, 0, 0, 10, 10, 10, tag=1)
    gmsh.model.occ.synchronize()

    # メッシュ生成
    gmsh.model.mesh.generate(3)  # 3は3次元メッシュを意味します

    # メッシュのエクスポート
    # Gmshのバージョンが4以降の場合、フォーマットバージョンを指定してエクスポートする必要があります。
    gmsh.option.setNumber("Mesh.MshFileVersion", 1.0)
    gmsh.write(output_filename)

    gmsh.finalize()

if __name__ == "__main__":
    create_3d_mesh("output_mesh.msh")
    # mesh = meshio.read('output_mesh.msh')   
    # # .vtkファイルに変換する
    # meshio.write("output_mesh.vtk", mesh, file_format="vtk")
    

    