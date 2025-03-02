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

def visualize_model(filename="model.msh"):
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

def create_3d_magnet_iron_mesh(filename="magnet_iron_model.msh"):
    gmsh.initialize()
    gmsh.model.add("MagnetIronModel")

    # 1x1x1の空間を定義
    lc = 0.1
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)
    p5 = gmsh.model.geo.addPoint(0, 0, 1, lc)
    p6 = gmsh.model.geo.addPoint(1, 0, 1, lc)
    p7 = gmsh.model.geo.addPoint(1, 1, 1, lc)
    p8 = gmsh.model.geo.addPoint(0, 1, 1, lc)

    # 磁石の大きさ、向き、重心座標を定義
    magnet_length = 0.3
    magnet_width = 0.1
    magnet_height = 0.05
    magnet_center = [0.5, 0.5, 0.4]  # 磁石の重心座標
    magnet_tilt = [math.radians(0), math.radians(0), math.radians(0)]  # 磁石の傾き（xyz成分）

    # 鉄の大きさ、向き、重心座標を定義
    iron_length = 0.3
    iron_width = 0.1
    iron_height = 0.05
    iron_center = [0.5, 0.5, 0.6]  # 鉄の重心座標
    iron_tilt = [math.radians(0), math.radians(0), math.radians(0)]  # 鉄の傾き（xyz成分）

    # 傾きを考慮して座標を計算する関数
    def calculate_tilted_point(x, y, z, center, tilt):
        # X軸回転
        y_rot = y * math.cos(tilt[0]) - z * math.sin(tilt[0])
        z_rot = y * math.sin(tilt[0]) + z * math.cos(tilt[0])
        y, z = y_rot, z_rot

        # Y軸回転
        x_rot = x * math.cos(tilt[1]) + z * math.sin(tilt[1])
        z_rot = -x * math.sin(tilt[1]) + z * math.cos(tilt[1])
        x, z = x_rot, z_rot

        # Z軸回転
        x_rot = x * math.cos(tilt[2]) - y * math.sin(tilt[2])
        y_rot = x * math.sin(tilt[2]) + y * math.cos(tilt[2])
        x, y = x_rot, y_rot

        return x, y, z

    # 磁石の頂点座標を計算
    pm1 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] - magnet_length / 2, magnet_center[1] - magnet_width / 2, magnet_center[2] - magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm2 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] + magnet_length / 2, magnet_center[1] - magnet_width / 2, magnet_center[2] - magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm3 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] + magnet_length / 2, magnet_center[1] + magnet_width / 2, magnet_center[2] - magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm4 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] - magnet_length / 2, magnet_center[1] + magnet_width / 2, magnet_center[2] - magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm5 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] - magnet_length / 2, magnet_center[1] - magnet_width / 2, magnet_center[2] + magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm6 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] + magnet_length / 2, magnet_center[1] - magnet_width / 2, magnet_center[2] + magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm7 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] + magnet_length / 2, magnet_center[1] + magnet_width / 2, magnet_center[2] + magnet_height / 2, magnet_center, magnet_tilt), lc)
    pm8 = gmsh.model.geo.addPoint(*calculate_tilted_point(magnet_center[0] - magnet_length / 2, magnet_center[1] + magnet_width / 2, magnet_center[2] + magnet_height / 2, magnet_center, magnet_tilt), lc)

    # 鉄の頂点座標を計算
    pi1 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] - iron_length / 2, iron_center[1] - iron_width / 2, iron_center[2] - iron_height / 2, iron_center, iron_tilt), lc)
    pi2 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] + iron_length / 2, iron_center[1] - iron_width / 2, iron_center[2] - iron_height / 2, iron_center, iron_tilt), lc)
    pi3 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] + iron_length / 2, iron_center[1] + iron_width / 2, iron_center[2] - iron_height / 2, iron_center, iron_tilt), lc)
    pi4 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] - iron_length / 2, iron_center[1] + iron_width / 2, iron_center[2] - iron_height / 2, iron_center, iron_tilt), lc)
    pi5 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] - iron_length / 2, iron_center[1] - iron_width / 2, iron_center[2] + iron_height / 2, iron_center, iron_tilt), lc)
    pi6 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] + iron_length / 2, iron_center[1] - iron_width / 2, iron_center[2] + iron_height / 2, iron_center, iron_tilt), lc)
    pi7 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] + iron_length / 2, iron_center[1] + iron_width / 2, iron_center[2] + iron_height / 2, iron_center, iron_tilt), lc)
    pi8 = gmsh.model.geo.addPoint(*calculate_tilted_point(iron_center[0] - iron_length / 2, iron_center[1] + iron_width / 2, iron_center[2] + iron_height / 2, iron_center, iron_tilt), lc)

    # 線を定義
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5 = gmsh.model.geo.addLine(p1, p5)
    l6 = gmsh.model.geo.addLine(p2, p6)
    l7 = gmsh.model.geo.addLine(p3, p7)
    l8 = gmsh.model.geo.addLine(p4, p8)
    l9 = gmsh.model.geo.addLine(p5, p6)
    l10 = gmsh.model.geo.addLine(p6, p7)
    l11 = gmsh.model.geo.addLine(p7, p8)
    l12 = gmsh.model.geo.addLine(p8, p5)

    # 磁石の線を定義
    lm1 = gmsh.model.geo.addLine(pm1, pm2)
    lm2 = gmsh.model.geo.addLine(pm2, pm3)
    lm3 = gmsh.model.geo.addLine(pm3, pm4)
    lm4 = gmsh.model.geo.addLine(pm4, pm1)
    lm5 = gmsh.model.geo.addLine(pm1, pm5)
    lm6 = gmsh.model.geo.addLine(pm2, pm6)
    lm7 = gmsh.model.geo.addLine(pm3, pm7)
    lm8 = gmsh.model.geo.addLine(pm4, pm8)
    lm9 = gmsh.model.geo.addLine(pm5, pm6)
    lm10 = gmsh.model.geo.addLine(pm6, pm7)
    lm11 = gmsh.model.geo.addLine(pm7, pm8)
    lm12 = gmsh.model.geo.addLine(pm8, pm5)

    # 鉄の線を定義
    li1 = gmsh.model.geo.addLine(pi1, pi2)
    li2 = gmsh.model.geo.addLine(pi2, pi3)
    li3 = gmsh.model.geo.addLine(pi3, pi4)
    li4 = gmsh.model.geo.addLine(pi4, pi1)
    li5 = gmsh.model.geo.addLine(pi1, pi5)
    li6 = gmsh.model.geo.addLine(pi2, pi6)
    li7 = gmsh.model.geo.addLine(pi3, pi7)
    li8 = gmsh.model.geo.addLine(pi4, pi8)
    li9 = gmsh.model.geo.addLine(pi5, pi6)
    li10 = gmsh.model.geo.addLine(pi6, pi7)
    li11 = gmsh.model.geo.addLine(pi7, pi8)
    li12 = gmsh.model.geo.addLine(pi8, pi5)

    # 面を定義
    ll1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    ll2 = gmsh.model.geo.addCurveLoop([l5, l9, -l6, -l1])
    ll3 = gmsh.model.geo.addCurveLoop([l6, l10, -l7, -l2])
    ll4 = gmsh.model.geo.addCurveLoop([l7, l11, -l8, -l3])
    ll5 = gmsh.model.geo.addCurveLoop([l8, l12, -l5, -l4])
    ll6 = gmsh.model.geo.addCurveLoop([l9, l10, l11, l12])

    # 磁石の面を定義
    llm1 = gmsh.model.geo.addCurveLoop([lm1, lm2, lm3, lm4])
    llm2 = gmsh.model.geo.addCurveLoop([lm5, lm9, -lm6, -lm1])
    llm3 = gmsh.model.geo.addCurveLoop([lm6, lm10, -lm7, -lm2])
    llm4 = gmsh.model.geo.addCurveLoop([lm7, lm11, -lm8, -lm3])
    llm5 = gmsh.model.geo.addCurveLoop([lm8, lm12, -lm5, -lm4])
    llm6 = gmsh.model.geo.addCurveLoop([lm9, lm10, lm11, lm12])

    # 鉄の面を定義
    lli1 = gmsh.model.geo.addCurveLoop([li1, li2, li3, li4])
    lli2 = gmsh.model.geo.addCurveLoop([li5, li9, -li6, -li1])
    lli3 = gmsh.model.geo.addCurveLoop([li6, li10, -li7, -li2])
    lli4 = gmsh.model.geo.addCurveLoop([li7, li11, -li8, -li3])
    lli5 = gmsh.model.geo.addCurveLoop([li8, li12, -li5, -li4])
    lli6 = gmsh.model.geo.addCurveLoop([li9, li10, li11, li12])

    s1 = gmsh.model.geo.addPlaneSurface([ll1])
    s2 = gmsh.model.geo.addPlaneSurface([ll2])
    s3 = gmsh.model.geo.addPlaneSurface([ll3])
    s4 = gmsh.model.geo.addPlaneSurface([ll4])
    s5 = gmsh.model.geo.addPlaneSurface([ll5])
    s6 = gmsh.model.geo.addPlaneSurface([ll6])

    # 磁石の面を定義
    sm1 = gmsh.model.geo.addPlaneSurface([llm1])
    sm2 = gmsh.model.geo.addPlaneSurface([llm2])
    sm3 = gmsh.model.geo.addPlaneSurface([llm3])
    sm4 = gmsh.model.geo.addPlaneSurface([llm4])
    sm5 = gmsh.model.geo.addPlaneSurface([llm5])
    sm6 = gmsh.model.geo.addPlaneSurface([llm6])

    # 鉄の面を定義
    si1 = gmsh.model.geo.addPlaneSurface([lli1])
    si2 = gmsh.model.geo.addPlaneSurface([lli2])
    si3 = gmsh.model.geo.addPlaneSurface([lli3])
    si4 = gmsh.model.geo.addPlaneSurface([lli4])
    si5 = gmsh.model.geo.addPlaneSurface([lli5])
    si6 = gmsh.model.geo.addPlaneSurface([lli6])

    # ボリュームを定義
    sl = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6])
    vol = gmsh.model.geo.addVolume([sl])

    # 磁石のボリュームを定義
    slm = gmsh.model.geo.addSurfaceLoop([sm1, sm2, sm3, sm4, sm5, sm6])
    volm = gmsh.model.geo.addVolume([slm])

    # 鉄のボリュームを定義
    sli = gmsh.model.geo.addSurfaceLoop([si1, si2, si3, si4, si5, si6])
    voli = gmsh.model.geo.addVolume([sli])

    # 物理グループを追加
    gmsh.model.addPhysicalGroup(3, [vol], name="Air")
    gmsh.model.addPhysicalGroup(3, [volm], name="Magnet")
    gmsh.model.addPhysicalGroup(3, [voli], name="Iron")

    # メッシュを生成
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # ファイルに書き出し
    gmsh.write(filename)
    gmsh.finalize()

if __name__ == "__main__":
    # 関数を呼び出して3Dモータモデルを作成
    # create_3d_motor_model(filename="model.msh")
    # visualize_model(filename="model.msh")
    create_3d_magnet_iron_mesh(filename="magnet_iron_model.msh")
    visualize_model(filename="magnet_iron_model.msh")