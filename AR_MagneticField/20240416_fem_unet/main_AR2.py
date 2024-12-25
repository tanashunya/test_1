import numpy as np
import cv2
from cv2 import aruco
import math
import gmsh
import Carray
from paraview.simple import *
import paraviewscript5 as para
import os
import time
from datetime import datetime
import config
# import MyFunctions2
import MyFunctions3 as myfunc3


class Marker:

    def __init__(self, id, corners, num_mag):
        self.id = id
        self.corners = corners
        self.vec = self.corners[2] - self.corners[1]
        self.thetarad = np.arctan2(self.vec[0], self.vec[1])
        self.thetadeg = math.degrees(self.thetarad)
        self.habaplus = self.corners[1][0] - self.corners[0][0]
        self.takasaplus = self.corners[1][1] - self.corners[0][1]
        self.vec_0to1 = self.corners[1] - self.corners[0]
        self.uppercenter = (self.corners[3] + self.corners[2]) / 2
        self.lowercenter = (self.corners[0] + self.corners[1]) / 2
        # print(self.uppercenter)

        self.virtual_corners = corners.copy()
        if config.VIRTUAL:
            self.virtual_corners[3][0] -= self.habaplus
            self.virtual_corners[3][1] -= self.takasaplus
            self.virtual_corners[2][0] += self.habaplus
            self.virtual_corners[2][1] += self.takasaplus
            self.virtual_corners[1][0] += self.habaplus
            self.virtual_corners[1][1] += self.takasaplus
            self.virtual_corners[0][0] -= self.habaplus
            self.virtual_corners[0][1] -= self.takasaplus
        
        # 1:空気 2:鉄 3:コイル 4:磁石
        if id == 0:
            self.material = 4 + num_mag     # 0 -> 磁石(4)
        else:
            self.material = 2     # それ以外： 鉄(2)
        # print(f"Marker ID: {self.id} theta: {self.thetadeg} material: {self.material}", self.id, "structed")

    # def __del__(self):
    #     print("Marker ID:", self.id, "deleted")
    def print_info(self):
        print(f"Marker ID: {self.id} theta: {self.thetarad} material: {self.material} ")

    def print_corners(self):
        print("Marker ID:", self.id)
        for i, corner in enumerate(self.corners):
            print("Corner", i+1, ":", corner)

    def set_loop_value(self, value):
        self.loop_value = value

    def set_plane_value(self, value):
        self.plane_value = value

    # def set_within(self):


def make_mesh(width, height,markers, lc):
    # メッシュ生成開始
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    mm = 1.0e-3

    # print(' ')
    # for marker in markers:
    #     marker.print_info()

    minus_loops = []

    # 空気以外の領域
    for marker in markers:
        if config.VIRTUAL:
            p1 = gmsh.model.geo.addPoint(marker.virtual_corners[3][0], height - (marker.virtual_corners[3][1]), 0.00, lc)
            p2 = gmsh.model.geo.addPoint(marker.virtual_corners[2][0], height - (marker.virtual_corners[2][1]), 0.00, lc)
            p3 = gmsh.model.geo.addPoint(marker.virtual_corners[1][0], height - (marker.virtual_corners[1][1]), 0.00, lc)
            p4 = gmsh.model.geo.addPoint(marker.virtual_corners[0][0], height - (marker.virtual_corners[0][1]), 0.00, lc)
        else:
            p1 = gmsh.model.geo.addPoint(marker.corners[3][0], height - (marker.corners[3][1]), 0.00, lc)
            p2 = gmsh.model.geo.addPoint(marker.corners[2][0], height - (marker.corners[2][1]), 0.00, lc)
            p3 = gmsh.model.geo.addPoint(marker.corners[1][0], height - (marker.corners[1][1]), 0.00, lc)
            p4 = gmsh.model.geo.addPoint(marker.corners[0][0], height - (marker.corners[0][1]), 0.00, lc)

        l1 = gmsh.model.geo.addLine(p1, p2)
        l2 = gmsh.model.geo.addLine(p2, p3)
        l3 = gmsh.model.geo.addLine(p3, p4)
        l4 = gmsh.model.geo.addLine(p4, p1)

        a_loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        marker.set_loop_value(gmsh.model.geo.addCurveLoop([l1, l2, l3, l4]))
        minus_loops.append(a_loop)
        
        marker.set_plane_value(gmsh.model.geo.addPlaneSurface([a_loop]))
    
    # 空気領域
    p1 = gmsh.model.geo.addPoint(0.00, 0.00, 0.00, lc)
    p2 = gmsh.model.geo.addPoint(width, 0.00, 0.00, lc)
    p3 = gmsh.model.geo.addPoint(width, height, 0.00, lc)
    p4 = gmsh.model.geo.addPoint(0.00, height, 0.00, lc)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    loop1 = []
    loop1.append(gmsh.model.geo.addCurveLoop([l1, l2, l3, l4]))
    # loop1のあとにminus_loopsをappend→plane1にいれる
    loop1.extend(minus_loops)
    # print(loop1)
    plane1 = gmsh.model.geo.addPlaneSurface(loop1)
    # まじない
    gmsh.model.geo.synchronize()

    pl2 = []
    pl3 = []
    pl4 = []
    pl5 = []
    pl6 = []
    pl7 = []
    pl8 = []
    pl9 = []
    pl10 = []
    pl11 = []
    pl12 = []
    pl13 = []

    for marker in markers:
        if marker.material == 2:
            pl2.append(marker.plane_value)
        elif marker.material == 3:
            pl3.append(marker.plane_value)
        elif marker.material == 4:
            pl4.append(marker.plane_value)
        elif marker.material == 5:
            pl5.append(marker.plane_value)
        elif marker.material == 6:
            pl6.append(marker.plane_value)
        elif marker.material == 7:
            pl7.append(marker.plane_value)
        elif marker.material == 8:
            pl8.append(marker.plane_value)
        elif marker.material == 9:
            pl9.append(marker.plane_value)
        elif marker.material == 10:
            pl10.append(marker.plane_value)
        elif marker.material == 11:
            pl11.append(marker.plane_value)
        elif marker.material == 12:
            pl12.append(marker.plane_value)
        elif marker.material == 13:
            pl13.append(marker.plane_value)
        
    group1  = gmsh.model.addPhysicalGroup(2, [plane1], 1) # 空気
    group2  = gmsh.model.addPhysicalGroup(2, pl2,      2) # 鉄
    group3  = gmsh.model.addPhysicalGroup(2, pl3,      3) # コイル
    group4  = gmsh.model.addPhysicalGroup(2, pl4,      4) # 磁石1
    group5  = gmsh.model.addPhysicalGroup(2, pl5,      5) # 磁石2
    group6  = gmsh.model.addPhysicalGroup(2, pl6,      6) # 磁石3
    group7  = gmsh.model.addPhysicalGroup(2, pl7,      7) # 磁石4
    group8  = gmsh.model.addPhysicalGroup(2, pl8,      8) # 磁石5
    group9  = gmsh.model.addPhysicalGroup(2, pl9,      9) # 磁石6
    group10 = gmsh.model.addPhysicalGroup(2, pl10,     10) # 磁石7
    group11 = gmsh.model.addPhysicalGroup(2, pl11,     11) # 磁石8
    group12 = gmsh.model.addPhysicalGroup(2, pl12,     12) # 磁石9
    group13 = gmsh.model.addPhysicalGroup(2, pl13,     13) # 磁石10

    # 材料の名前つけ
    gmsh.model.setPhysicalName(2, group1, "Air")
    gmsh.model.setPhysicalName(2, group2, "Iron")
    gmsh.model.setPhysicalName(2, group3, "Coil")
    gmsh.model.setPhysicalName(2, group4, "Magnet1")
    gmsh.model.setPhysicalName(2, group5, "Magnet2")
    gmsh.model.setPhysicalName(2, group6, "Magnet3")
    gmsh.model.setPhysicalName(2, group7, "Magnet4")
    gmsh.model.setPhysicalName(2, group8, "Magnet5")
    gmsh.model.setPhysicalName(2, group9, "Magnet6")
    gmsh.model.setPhysicalName(2, group10, "Magnet7")
    gmsh.model.setPhysicalName(2, group11, "Magnet8")
    gmsh.model.setPhysicalName(2, group12, "Magnet9")
    gmsh.model.setPhysicalName(2, group13, "Magnet10")
    
    # メッシュ生成（2次元）
    if gmsh.model.mesh.generate(2)==False:
        # print("メッシュ作れなかったよ。")
        # gmsh.model.removeEntities()
        # gmsh.model.removePhysicalGroups()
        # gmsh.initialize()
        # gmsh.write("gomi.msh1")
        # gmsh.finalize()
        return False

    gmsh.write("variablemesh2.msh1")
    # gmsh.write("meshdata" + datetime.now().strftime("%Y%m%d_%H%M%S") + ".msh1")

    gmsh.finalize()

def make_input_data(markers, path):
    # 学習データ用の白い背景画像を用意
    InputDataImage = np.zeros((config.HEIGHT, config.WIDTH, 3),dtype=np.uint8 ) + 255
    b1, g1, r1 = 255, 0, 0   # 青
    b2, g2, r2 = 0, 0, 255   # 赤
    b, g, r    = 0, 255, 0   # 緑
    for marker in markers:
        # print(marker.id)
        if marker.id == id_magnet:
            # points1 = np.array([marker.virtual_corners[3],marker.corners[2]-(marker.vec_0to1/2),marker.corners[1]-(marker.vec_0to1/2),marker.virtual_corners[0]])
            # points2 = np.array([marker.corners[2]-(marker.vec_0to1/2),marker.virtual_corners[2],marker.virtual_corners[1],marker.corners[1]-(marker.vec_0to1/2)])
            points1 = np.array([[marker.virtual_corners[0], marker.lowercenter, marker.uppercenter, marker.virtual_corners[3]]])
            points1 = points1.astype(int)
            points2 = np.array([[marker.lowercenter, marker.virtual_corners[1], marker.virtual_corners[2], marker.uppercenter]])
            points2 = points2.astype(int)
            # print(points1)
            
            InputDataImage = cv2.fillPoly(InputDataImage, points1, color=(b1, g1, r1))
            InputDataImage = cv2.fillPoly(InputDataImage, points2, color=(b2, g2, r2))
        else:
            points = np.array([[marker.virtual_corners[0],marker.virtual_corners[1],marker.virtual_corners[2],marker.virtual_corners[3]]])
            points = points.astype(int)
            InputDataImage = cv2.fillPoly(InputDataImage, points, color=(b, g, r))
    cv2.imwrite(path, InputDataImage)

# 縦幅と横幅（解像度）
imgsize_x = config.WIDTH
imgsize_y = config.HEIGHT
height = float(imgsize_y)
width = float(imgsize_x)

# ウィンドウの設定
window_name = 'press q to exit'
x_position = 640
y_position = 0
# ウィンドウを作成
cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
cv2.moveWindow(window_name, x_position, y_position)
# cv2.resizeWindow(window_name, 1280, 960)

# カメラの設定
cap = cv2.VideoCapture(0)  # 0はデフォルトのウェブカメラ

# ARマーカの設定
aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_50)
aruco_params = aruco.DetectorParameters_create()
id_magnet = 0

# mesh size
lc = config.MESH_SIZE

# 現在のディレクトリを取得
current_dir = os.getcwd()

while True:
    loopstarttime = time.perf_counter()
    markers = []
    ret, frame = cap.read()  # フレームのキャプチャ
    if config.WIDTH!=640 or config.HEIGHT!=480:
        frame = frame[0:config.WIDTH, 0:config.HEIGHT]
    if config.ROTATE_180:
        frame = cv2.flip(frame, -1) # 180度回す
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)  # グレースケールへの変換
    corners, ids, rejected = aruco.detectMarkers(gray, aruco_dict, parameters=aruco_params)  # ARマーカの検出

    # マーカが検出された場合、インスタンスを生成
    if ids is not None:
        num_mag = 0
        num_iron = 0
        for i in range(len(ids)):
            marker_id = ids[i][0]
            marker_corners = corners[i][0]
            # marker.print_corners()
            m = Marker(marker_id, marker_corners, num_mag)
            markers.append(m)
            if myfunc3.is_sticking_out(m.virtual_corners):
                markers.pop(-1)
                print("はみ出しているオブジェクトを削除しました。")
            elif marker_id == id_magnet:
                num_mag += 1
            elif marker_id == 2:
                num_iron += 1

        # print(' ')
        # for marker in markers:
        #     marker.print_info()
        
        if num_mag==0:
            print("磁石がありません")
        elif myfunc3.is_possible_to_make_mesh(markers) is False and config.VIRTUAL:
            print("メッシュを生成できませんでした")
        else:
            path_of_output_data = ""
            if config.DROW_OBJECT_CONTOUR_MODE:
                frame = myfunc3.draw_obj_contour(frame, markers)
            if config.MAKE_DATASET:
                i = config.INITIAL_INDEX
                while(os.path.exists(config.DIR_OF_INPUT_DATA+f"{i}.png")):
                    i += 1
                print(f"data_index = {i}")
                path_of_input_data = config.DIR_OF_INPUT_DATA+f"{i}.png"
                path_of_output_data = config.DIR_OF_OUTPUT_DATA+f"{i}.csv"
                make_input_data(markers, path_of_input_data)
            make_mesh(width, height, markers, lc)
            femstart = time.perf_counter()
            # FEM ----------------------------------------------------------------
            # for marker in markers:
            #     if marker.material==
            
            theta_of_magnets = np.array([])
            for marker in markers:
                if marker.id == id_magnet:
                    theta_of_magnets = np.append(theta_of_magnets, marker.thetarad)
            # print(theta_of_magnets)
            theta = 0.0
            my_carray = Carray.Carray(theta, 
                                      len(theta_of_magnets), 
                                      theta_of_magnets, 
                                      int(width), int(height), 
                                      config.MAKE_DATASET, i, 
                                      path_of_output_data)
            my_carray.fem()
            # print(time.perf_counter() - femstart)
            # ベクトル or コンター描画--------------------------------------------------------
            if config.DRAW_VECTOR:
                # frame = Myfunctions.draw_vectors(frame, imgsize_y) # Legacy ver.
                para.makevecpng(current_dir)
                overlay_image = cv2.imread("paraview_img.png")
                frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)
            elif config.DRAW_CONTOUR:
                para.makeconpng(current_dir)
                overlay_image = cv2.imread("contour.png")
                frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)


    cv2.imshow(window_name, frame)


    if cv2.waitKey(1) == ord('q'):  # qを押すと終了
        if config.FRAME_SAVE_MODE:
            cv2.imwrite("frame.png", frame)
            print("frame.png saved successfully!")
        break
    
    print(time.perf_counter() - loopstarttime)


cap.release()  # カメラの解放
cv2.destroyAllWindows()  # ウィンドウを閉じる
