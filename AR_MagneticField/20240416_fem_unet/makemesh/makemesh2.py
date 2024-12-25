import gmsh

def create_cube_mesh():
    gmsh.initialize()
    gmsh.model.add("Cube")
    # 六面体要素でのメッシュ生成を指示
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    # 立方体の頂点を定義
    lc = 0.5  # メッシュの局所的な特徴サイズ
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(1, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(1, 1, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, 1, 0, lc)
    p5 = gmsh.model.geo.addPoint(0, 0, 1, lc)
    p6 = gmsh.model.geo.addPoint(1, 0, 1, lc)
    p7 = gmsh.model.geo.addPoint(1, 1, 1, lc)
    p8 = gmsh.model.geo.addPoint(0, 1, 1, lc)

    # 立方体の辺を定義
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)
    l9 = gmsh.model.geo.addLine(p1, p5)
    l10 = gmsh.model.geo.addLine(p2, p6)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p4, p8)

    # 立方体の面を形成するループを定義
    cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
    cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
    cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
    cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])

    # 立方体の面を定義
    s1 = gmsh.model.geo.addPlaneSurface([cl1])
    s2 = gmsh.model.geo.addPlaneSurface([cl2])
    s3 = gmsh.model.geo.addPlaneSurface([cl3])
    s4 = gmsh.model.geo.addPlaneSurface([cl4])
    s5 = gmsh.model.geo.addPlaneSurface([cl5])
    s6 = gmsh.model.geo.addPlaneSurface([cl6])

    # サーフェスループを定義
    sl = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6])

    # 立方体の体積を定義
    vol = gmsh.model.geo.addVolume([sl])

    gmsh.model.geo.synchronize()

    # # Transfiniteアルゴリズムを使用して六面体メッシュを生成
    # gmsh.model.mesh.setTransfiniteSurface(s1)
    # gmsh.model.mesh.setTransfiniteSurface(s2)
    # gmsh.model.mesh.setTransfiniteSurface(s3)
    # gmsh.model.mesh.setTransfiniteSurface(s4)
    # gmsh.model.mesh.setTransfiniteSurface(s5)
    # gmsh.model.mesh.setTransfiniteSurface(s6)

    # gmsh.model.mesh.setTransfiniteVolume(vol)

    # メッシュ生成
    gmsh.model.mesh.generate(3)

    # ファイルに保存
    gmsh.write("cube3.msh1")

    gmsh.finalize()

if __name__ == "__main__":
    create_cube_mesh()