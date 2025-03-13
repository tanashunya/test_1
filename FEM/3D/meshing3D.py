"""
このスクリプトは、Gmshライブラリを使用してローターとステーターのジオメトリのメッシュモデルを生成します。
メッシュはディスクに書き込まれ、その後、mapdl-archiveを使用して読み取られ、アーカイブファイルに変換され、
さらにFEM解析を行います。このモジュールは、オプションの可視化のためにpyvistaを使用します。
単一の周期セクターを生成します。
必要条件
------------
.. code::
   pip install gmsh mapdl-archive pyvista numpy
"""
import math
import time
import sys

import gmsh
import mapdl_archive
import numpy as np
import pyvista as pv
import meshio

def create_3d_motor_model(filename="model.msh"):
    # 幾何学的パラメータを定義
    r1 = 0.057  # 内部半径
    r2 = 0.067  # ローター溝の底から中心まで
    r3 = 0.08  # ローター溝の上部から中心まで
    r4 = 0.0814  # ステーター溝の上部から中心まで
    r5 = 0.0944  # ローター溝の底から中心まで
    r6 = 0.1044  # 外部半径（直径の半分）

    N_ROTOR = 400
    N_STATOR = 440

    PT_ROTOR = 0.5  # ローター歯の充填率（0から1まで）
    PT_STATOR = 0.5  # ステーター歯の充填率（0から1まで）

    # 状態
    AR = 2  # 機械的度数でのローター位置

    N_RAD_DENSITY = 4
    N_TAN_DENSITY = 10

    def polar_to_cartesian(radius, angle, center):
        x = center[0] + radius * math.cos(angle)
        y = center[1] + radius * math.sin(angle)
        return (x, y)

    # Gmsh APIを開始
    gmsh.initialize()
    gmsh.model.add("Model")

    # gmshモデルにポイントを追加
    center_xy = (0.0, 0.0)
    p_center = gmsh.model.geo.addPoint(*center_xy, 0.0)

    def gen_arcs(n_arc, radius, ang_start=0, fill_factor=0.5):
        """すべての内側の弧を生成します。"""
        arc_len = math.pi / (n_arc / 2)
        points = []
        angle = ang_start
        for ii in range(n_arc):
            x, y = polar_to_cartesian(radius, angle, center_xy)
            if ii % 2:
                angle += 2 * arc_len * fill_factor
            else:
                angle += 2 * arc_len * (1 - fill_factor)
            points.append(gmsh.model.geo.addPoint(x, y, 0.0))

        # 弧を作成
        arcs = []
        for ii in range(n_arc - 1):
            arcs.append(gmsh.model.geo.addCircleArc(points[ii], p_center, points[ii + 1]))

        arcs.append(gmsh.model.geo.addCircleArc(points[-1], p_center, points[0]))

        # トランスフィニットラインを設定
        for ii in range(n_arc):
            gmsh.model.geo.mesh.setTransfiniteCurve(arcs[ii], N_TAN_DENSITY)

        return points, arcs

    def gen_cyc(n_arc, r_inner, r_middle, ang_start=0.0, fill_factor=0.5):
        ang_start_rad = ang_start * math.pi / 180
        points_inner, arcs_inner = gen_arcs(n_arc, r_inner, ang_start_rad, fill_factor)
        loop_inner = gmsh.model.geo.addCurveLoop(arcs_inner)

        points_middle, arcs_middle = gen_arcs(n_arc, r_middle, ang_start_rad, fill_factor)
        loop_middle = gmsh.model.geo.addCurveLoop(arcs_middle)

        # 内側から中間への接続ラインを生成
        in_to_mid_con_lines = []
        for ii in range(n_arc):
            in_to_mid_con_lines.append(
                gmsh.model.geo.addLine(points_inner[ii], points_middle[ii])
            )
            gmsh.model.geo.mesh.setTransfiniteCurve(in_to_mid_con_lines[-1], N_RAD_DENSITY)

        # 磁石と鉄心を表す個々の「ループ」を生成
        plane_inner = []
        # plane_outer = []
        # for ii in range(n_arc):
        for ii in range(1):  # 単一の弧のオーバーライド
            # 「内側のコア」を生成、これらはすべて鉄です
            line0 = in_to_mid_con_lines[ii]
            # 最後のラインをラップアラウンド
            if ii == n_arc - 1:
                line1 = in_to_mid_con_lines[0]
            else:
                line1 = in_to_mid_con_lines[ii + 1]
            loop = gmsh.model.geo.addCurveLoop(
                [arcs_inner[ii], line0, arcs_middle[ii], line1], reorient=True
            )
            plane = gmsh.model.geo.addPlaneSurface([loop])
            plane_inner.append(plane)
            gmsh.model.geo.mesh.setTransfiniteSurface(plane)

        return plane_inner

    ###############################################################################
    # 2Dジオメトリを生成
    inner_plane = gen_cyc(4, 1, 2)
    gmsh.model.geo.synchronize()

    # 四角形要素を生成するオプションを設定
    gmsh.option.setNumber("Mesh.RecombineAll", True)

    # 2Dメッシュを生成
    # gmsh.model.mesh.generate(2)

    ###############################################################################
    # 3Dジオメトリを生成
    # 押し出しベクトルを定義 [dx, dy, dz]
    extrude_vector = [0, 0, 0.5]  # 押し出し長さを0.1に置き換え

    pgroup = 2
    base_surf_pg = gmsh.model.addPhysicalGroup(
        pgroup, inner_plane, tag=100, name="lower_surface"
    )

    # メッシュを押し出し
    subdivision = [3]
    dimTags = gmsh.model.getEntities(2)
    extrusion = gmsh.model.geo.extrude(
        [(pgroup, inner_plane[0])], 0, 0, 1, subdivision, recombine=True
    )

    ###############################################################################
    # 3Dメッシュ
    gmsh.model.geo.synchronize()
    volume = gmsh.model.addPhysicalGroup(3, [extrusion[1][1]], name="volume")

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.refine()

    ###############################################################################
    # ファイルに書き出し
    filename = "model.msh"
    gmsh.write(filename)

def visualize_mesh_by_pyvista(filename="model.msh"):
    ###############################################################################
    # pyvistaを使用して読み込み、mapdl-archiveを使用してアーカイブファイルに変換
    grid = pv.read(filename)
    # grid.plot(color='w', show_edges=True, line_width=1.5)

    mapdl_archive.save_as_archive("archive.cdb", grid)

    ###############################################################################
    # 可視化する
    plotter = pv.Plotter()
    plotter.add_axes()  # xyz軸を表示させる
    plotter.add_mesh(grid, color='w', show_edges=True, line_width=1.5, opacity=0.5)  # 空気領域を半透明に設定
    plotter.show()

class MaterialBlock:
    """直方体を表すクラス"""
    def __init__(self, center, size, rotation=(0, 0, 0), name="Block", lc=0.1):
        self.center = center
        self.size = size
        self.rotation = [math.radians(r) if isinstance(r, (int, float)) else r for r in rotation]
        self.name = name
        self.lc = lc

    def create(self):
        # 中心を基準に直方体を作成
        x, y, z = [c - s / 2 for c, s in zip(self.center, self.size)]
        tag = gmsh.model.occ.addBox(x, y, z, *self.size)
        
        # 回転を適用
        if any(self.rotation):
            gmsh.model.occ.rotate([(3, tag)], self.center[0], self.center[1], self.center[2], 1, 0, 0, self.rotation[0])
            gmsh.model.occ.rotate([(3, tag)], self.center[0], self.center[1], self.center[2], 0, 1, 0, self.rotation[1])
            gmsh.model.occ.rotate([(3, tag)], self.center[0], self.center[1], self.center[2], 0, 0, 1, self.rotation[2])
        
        return tag

def check_intersection(air_tag, material_tags):
    """材料同士の交差をチェックし、交差があればエラーメッセージを出力して終了"""
    for i, tag1 in enumerate(material_tags):
        for tag2 in material_tags[i+1:]:
            gmsh.model.occ.fragment([(3, tag1)], [(3, tag2)])
            gmsh.model.occ.synchronize()
            entities = gmsh.model.getEntities(3)
            if len(entities) > len(material_tags) + 1:  # 交差がある場合
                print(f"Error: Materials {tag1} and {tag2} are intersecting.")
                gmsh.finalize()
                sys.exit(1)
    gmsh.model.occ.synchronize()

def create_3d_magnet_iron_mesh(filename):
    """磁石と鉄心のメッシュモデルを生成"""
    gmsh.initialize()
    gmsh.model.add("MagnetIronModel")

    # 空気領域の直方体
    air = MaterialBlock(
        center=[0.5, 0.5, 0.5],
        size=[1.0, 1.0, 1.0],
        rotation=[0, 0, 0],
        name="Air",
        lc=0.1
    )
    air_tag = air.create()

    # 磁石の直方体
    magnet = MaterialBlock(
        center=[0.5, 0.5, 0.4],
        size=[0.3, 0.1, 0.05],
        rotation=[0, 0, 0],
        name="Magnet",
        lc=0.1
    )
    magnet_tag = magnet.create()

    # 鉄心の直方体
    iron = MaterialBlock(
        center=[0.5, 0.5, 0.6],
        size=[0.3, 0.1, 0.05],
        rotation=[0, 0, 0],
        name="Iron",
        lc=0.1
    )
    iron_tag = iron.create()

    # 材料同士の交差をチェック（空気は除外）
    check_intersection(air_tag, [magnet_tag, iron_tag])

    # 材料と空気領域のブール演算 (差分)
    gmsh.model.occ.cut([(3, air_tag)], [(3, magnet_tag), (3, iron_tag)])
    gmsh.model.occ.synchronize()

    # 空気領域の外部面に物理グループを追加
    air_faces = gmsh.model.getBoundary([(3, air_tag)], oriented=False)
    gmsh.model.addPhysicalGroup(2, [tag for dim, tag in air_faces if dim == 2], name="Boundary")

    # 四面体メッシュを生成するための設定
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # 1はDelaunay法を指定

    # メッシュを生成
    gmsh.model.mesh.generate(3)

    # ファイルに書き出し
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    # 関数を呼び出して3Dモータモデルを作成
    # create_3d_motor_model(filename="model.msh")
    # visualize_mesh_by_pyvista(filename="model.msh")
    create_3d_magnet_iron_mesh(filename="3d_magnet_iron_model.msh")
    create_3d_magnet_iron_mesh(filename="3d_magnet_iron_model.msh1")
    visualize_mesh_by_pyvista(filename="3d_magnet_iron_model.msh")