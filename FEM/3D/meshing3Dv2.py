import gmsh
from meshing3D import visualize_mesh_by_pyvista

# Gmshの初期化
gmsh.initialize()

# 新しいモデルを作成
gmsh.model.add("cube")

# 1x1x1の立方体を定義
lc = 0.1  # メッシュサイズ
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)

# 形状を同期
gmsh.model.occ.synchronize()

# メッシュサイズを設定
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc)

# 6面体メッシュを生成するための設定
gmsh.option.setNumber("Mesh.Algorithm3D", 8)  # 8はHXT
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # 1はHexahedral

# メッシュを生成
gmsh.model.mesh.generate(3)

# メッシュをファイルに保存
gmsh.write("cube.msh")
gmsh.write("cube.msh1")

# Gmshを終了
gmsh.finalize()

visualize_mesh_by_pyvista(filename="cube.msh")
