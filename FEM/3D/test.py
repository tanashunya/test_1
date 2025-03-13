import gmsh
import numpy as np
import pyvista as pv
import mapdl_archive


class MaterialBlock:
    def __init__(self, x, y, z, dx, dy, dz, rotation=0):
        self.x = x
        self.y = y
        self.z = z
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rotation = rotation

    def create(self):
        tag = gmsh.model.occ.addBox(self.x, self.y, self.z, self.dx, self.dy, self.dz)
        if self.rotation != 0:
            gmsh.model.occ.rotate([(3, tag)], self.x, self.y, self.z, 0, 0, 1, self.rotation)
        return tag

def create_mesh(meshname, materials):
    gmsh.initialize()
    gmsh.model.add("3D_Mesh")

    # 空気領域の作成 (1x1x1)
    air_tag = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)

    # 各材料の配置
    material_tags = [material.create() for material in materials]

    # 材料と空気領域のブール演算 (差分)
    for tag in material_tags:
        gmsh.model.occ.cut([(3, air_tag)], [(3, tag)])

    gmsh.model.occ.synchronize()

    # メッシュの生成
    gmsh.model.mesh.generate(3)

    # メッシュファイルの書き出し
    gmsh.write(meshname)

    gmsh.finalize()

def visualize_mesh(meshname):
    ###############################################################################
    # pyvistaを使用して読み込み、mapdl-archiveを使用してアーカイブファイルに変換
    grid = pv.read(meshname)
    # grid.plot(color='w', show_edges=True, line_width=1.5)

    mapdl_archive.save_as_archive("archive.cdb", grid)

    ###############################################################################
    # 可視化する
    plotter = pv.Plotter()
    plotter.add_axes()  # xyz軸を表示させる
    plotter.add_mesh(grid, color='w', show_edges=True, line_width=1.5, opacity=0.5)  # 空気領域を半透明に設定
    plotter.show()


if __name__ == "__main__":
    materials = [
        MaterialBlock(0.2, 0.2, 0.2, 0.4, 0.4, 0.4, rotation=np.pi / 6),
        MaterialBlock(0.5, 0.5, 0.5, 0.3, 0.3, 0.3, rotation=np.pi / 4)
    ]
    create_mesh(meshname="mesh.msh", materials=materials)
    create_mesh(meshname="mesh.msh1", materials=materials)
    visualize_mesh(meshname="mesh.msh")
