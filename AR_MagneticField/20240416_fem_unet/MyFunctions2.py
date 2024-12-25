import gmsh
import csv
import math
import cv2
import numpy as np
import time
# import mygmshfunc
#############################################################################################
def generate_mesh(width, height, corners):
    # メッシュ生成開始
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    mm = 1.0e-3

    # メッシュの間隔の定義
    lc = 20.0
    lc2 = 20.0
    
    # 点
    p1 = gmsh.model.geo.addPoint(0.00, 0.00, 0.00, lc)
    p2 = gmsh.model.geo.addPoint(width, 0.00, 0.00, lc)
    p3 = gmsh.model.geo.addPoint(width, height, 0.00, lc)
    p4 = gmsh.model.geo.addPoint(0.00, height, 0.00, lc)

    p5 = gmsh.model.geo.addPoint(corners[0][0][3][0], height - corners[0][0][3][1], 0.00, lc2)
    p6 = gmsh.model.geo.addPoint(corners[0][0][2][0], height - corners[0][0][2][1], 0.00, lc2)
    p7 = gmsh.model.geo.addPoint(corners[0][0][1][0], height - corners[0][0][1][1], 0.00, lc2)
    p8 = gmsh.model.geo.addPoint(corners[0][0][0][0], height - corners[0][0][0][1], 0.00, lc2)

    # 線
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)

    # 閉路
    loop1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    loop2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])

    # 閉曲面
    plane1 = gmsh.model.geo.addPlaneSurface([loop1, loop2])  # loop1からloop2を引いた面をplane1にした
    plane2 = gmsh.model.geo.addPlaneSurface([loop2])

    # シンクロナイズ
    gmsh.model.geo.synchronize()

    # グループ（次元数、[閉曲面]、材料番号）
    group1 = gmsh.model.addPhysicalGroup(2, [plane1], 1)  # 空気
    group2 = gmsh.model.addPhysicalGroup(2, [plane2], 4)  # 磁石

    # 材料の名前つけ
    gmsh.model.setPhysicalName(2, group1, "Air")
    gmsh.model.setPhysicalName(2, group2, "Magnet")

    # メッシュ生成（2次元）
    gmsh.model.mesh.generate(2)

    # メッシュデータを配列として取得
    nodeTags,coord,parametricCoord = gmsh.model.mesh.getNodes()
    elementTypes,elementTags,nodeTagsE = gmsh.model.mesh.getElements()
    # ノード数print
    num_nodes = len(nodeTags)
    # print(num_nodes)
    # エレメント数print
    num_elements = len(elementTags[1])
    # print(num_elements)
    # print(coord)
    # 二次元に変換
    # coord_2d = [coord[i:i+3] for i in range(0, num_elements, 3)]
    # coord_2d = coord_2d.T
    # print(coord_2d)

    gmsh.write("variablemesh.msh1")
    # gmsh.write("meshdata" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".msh1")

    gmsh.finalize()
#############################################################################################





#############################################################################################
def draw_vectors(frame, imgsize_y):
    readcsvstart = time.perf_counter()#####
    vectors = []
    length = []
    truelength = []

    with open('vectors.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            x_coordinate = float(row[0])
            y_coordinate = float(imgsize_y) - float(row[1])
            x_component = float(row[3]) * 6
            y_component = -float(row[4]) * 6
            truelength.append(math.sqrt(float(row[3]) ** 2 + float(row[4]) ** 2))
            length.append(math.sqrt(x_component ** 2 + y_component ** 2))
            vectors.append((x_coordinate, y_coordinate, x_component, y_component))

    # print(time.perf_counter() - readcsvstart)#####
    drawvecstart = time.perf_counter()#####
    colormap = cv2.COLORMAP_JET  # 使用するカラーマップを指定

    min_length = min(length)
    max_length = max(length)
    min_true_length = min(truelength)
    max_true_length = max(truelength)

    # frame = np.ones((480, 640, 3), dtype=np.uint8) * 255  # 描画用の白い背景画像を作成frame = np.ones((480, 640, 3), dtype=np.uint8) * 255  # 描画用の白い背景画像を作成

    for i, vector in enumerate(vectors):
        if(i%2==0):
            x_coordinate, y_coordinate, x_component, y_component = vector
            vector_length = math.sqrt(x_component ** 2 + y_component ** 2)

            # ベクトルの大きさを最小値から最大値に正規化
            normalized_length = (vector_length - min_length) / (max_length - min_length)

            # カラーマップから色を取得
            # grayscale_value = np.uint8(np.array([normalized_length]) * 255)

            # 正規化されたベクトルの大きさをグレースケールに変換
            grayscale_value = np.uint8(np.array([normalized_length]) * 255)

            # グレースケール画像にカラーマップを適用
            color = cv2.applyColorMap(grayscale_value, colormap)

            scaled_x_component = x_component / vector_length * 6
            scaled_y_component = y_component / vector_length * 6

            start_point = (int(x_coordinate), int(y_coordinate))
            end_point = (int(x_coordinate + scaled_x_component), int(y_coordinate + scaled_y_component))

            frame = cv2.arrowedLine(frame, start_point, end_point,(int(color[0,0,0]), int(color[0,0,1]),int(color[0,0,2])), thickness=1)

    # print(time.perf_counter()-drawvecstart)#####

    return frame
#############################################################################################


