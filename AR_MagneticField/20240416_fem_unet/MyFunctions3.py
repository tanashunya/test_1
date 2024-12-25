import config
import numpy as np
import cv2
import torch
import gmsh
import math
import random
import sys
import os
import Carray


class Random_Object:
    def __init__(self, size):
        x = np.random.randint(config.WIDTH)
        y = np.random.randint(config.HEIGHT)
        theta = random.uniform(0, 2*math.pi)
        # center = np.array([x, y])
        
        # 正方形の各頂点の相対座標を計算
        vertices_relative = np.array([
            [size / 2, size / 2],
            [-size / 2, size / 2],
            [-size / 2, -size / 2],
            [size / 2, -size / 2]
        ])
        # 回転行列を作成し、正方形の頂点を計算
        rotation_matrix = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ])
        vertices_rotated = np.dot(vertices_relative, rotation_matrix)
    
        # 中心座標を加えて絶対座標を計算
        self.coordinates = vertices_rotated + np.array([x,y])

        # self.coordinates = np.array([[x-size/2, y+size/2],
        #                         [x+size/2, y+size/2],
        #                         [x+size/2, y-size/2],
        #                         [x-size/2, y-size/2]])
        # print(self.coordinates)
        # # 回転行列
        # rotation_matrix = np.array([
        #     [np.cos(theta), np.sin(theta)],
        #     [np.sin(theta), np.cos(theta)]
        # ])
        # # 中心座標を原点に移動　回転　再び中心座標の位置に戻す
        # centered_coordinates = coordinates - center
        # self.coordinates = np.dot(centered_coordinates, rotation_matrix) + center
        # print(self.coordinates)

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

def calculate_square_vertices(center_x, center_y, side_length, angle_deg):
    # 角度をラジアンに変換
    angle_rad = np.radians(angle_deg)
    
    # 正方形の各頂点の相対座標を計算
    vertices_relative = np.array([
        [side_length / 2, side_length / 2],
        [-side_length / 2, side_length / 2],
        [-side_length / 2, -side_length / 2],
        [side_length / 2, -side_length / 2]
    ])
    
    # 回転行列を作成し、正方形の頂点を計算
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad), np.cos(angle_rad)]
    ])
    vertices_rotated = np.dot(vertices_relative, rotation_matrix)
    
    # 中心座標を加えて絶対座標を計算
    vertices = vertices_rotated + np.array([center_x, center_y])
    
    return vertices

def is_possible_to_make_mesh(markers):
    # 衝突判定
    for marker in markers:
        if config.VIRTUAL:
            if is_sticking_out(marker.virtual_corners):
                # print("オブジェクトがフレームからはみ出ています")
                return False
            for m in markers:
                if marker.id != m.id:
                    qs1 = [(marker.virtual_corners[0][0], marker.virtual_corners[0][1]), (marker.virtual_corners[1][0], marker.virtual_corners[1][1]), (marker.virtual_corners[2][0], marker.virtual_corners[2][1]), (marker.virtual_corners[3][0], marker.virtual_corners[3][1])]
                    qs2 = [(m.virtual_corners[0][0], m.virtual_corners[0][1]), (m.virtual_corners[1][0], m.virtual_corners[1][1]), (m.virtual_corners[2][0], m.virtual_corners[2][1]), (m.virtual_corners[3][0], m.virtual_corners[3][1])]
                    # print(convex_polygons_intersection(qs1, qs2))
                    if convex_polygons_intersection(qs1, qs2)==1:
                        # print("オブジェクト同士が重なっています")
                        return False
    return True

def dot3(O, A, B):
    ox, oy = O; ax, ay = A; bx, by = B
    return (ax - ox) * (bx - ox) + (ay - oy) * (by - oy)
def cross3(O, A, B):
    ox, oy = O; ax, ay = A; bx, by = B
    return (ax - ox) * (by - oy) - (bx - ox) * (ay - oy)
def dist2(A, B):
    ax, ay = A; bx, by = B
    return (ax - bx) ** 2 + (ay - by) ** 2
def is_intersection(P0, P1, Q0, Q1):
    C0 = cross3(P0, P1, Q0)
    C1 = cross3(P0, P1, Q1)
    if C0 == C1 == 0:
        E0 = dot3(P0, P1, Q0)
        E1 = dot3(P0, P1, Q1)
        if not E0 < E1:
            E0, E1 = E1, E0
        return 0 <= E1 and E0 <= dist2(P0, P1)
    D0 = cross3(Q0, Q1, P0)
    D1 = cross3(Q0, Q1, P1)
    return C0 * C1 <= 0 and D0 * D1 <= 0
def line_cross_point(P0, P1, Q0, Q1):
    x0, y0 = P0; x1, y1 = P1
    x2, y2 = Q0; x3, y3 = Q1
    dx0 = x1 - x0; dy0 = y1 - y0
    dx1 = x3 - x2; dy1 = y3 - y2

    s = (y0-y2)*dx1 - (x0-x2)*dy1
    sm = dx0*dy1 - dy0*dx1
    return (x0 + s*dx0/sm, y0 + s*dy0/sm) if s != 0 else (x0, y0)
# ps, qs: a polygon (counter-clockwise)
def convex_polygons_intersection(ps, qs):
    pl = len(ps); ql = len(qs)
    i = j = 0
    while (i < pl or j < ql) and (i < 2*pl) and (j < 2*ql):
        px0, py0 = ps0 = ps[(i-1)%pl]; px1, py1 = ps1 = ps[i%pl]
        qx0, qy0 = qs0 = qs[(j-1)%ql]; qx1, qy1 = qs1 = qs[j%ql]

        if is_intersection(ps0, ps1, qs0, qs1):
            return 1

        ax = px1 - px0; ay = py1 - py0
        bx = qx1 - qx0; by = qy1 - qy0

        v = (ax*by - bx*ay)
        va = cross3(qs0, qs1, ps1)
        vb = cross3(ps0, ps1, qs1)

        if v == 0 and va < 0 and vb < 0:
            return 0
        if v == 0 and va == 0 and vb == 0:
            i += 1
        elif v >= 0:
            if vb > 0:
                i += 1
            else:
                j += 1
        else:
            if va > 0:
                j += 1
            else:
                i += 1
    return 0

def is_sticking_out(corners):
    # for corner in corners:
    #     print(corner)
    for corner in corners:
        if corner[0] <= 0 or config.WIDTH <=corner[0]: # x座標のチェック
            return True
        if corner[1] <= 0 or config.HEIGHT <=corner[1]: # x座標のチェック
            return True
    return False

# 四角形しかできない、だめじゃん
def draw_obj_contour(frame, markers):
    for marker in markers:
        if config.VIRTUAL:
            points = np.array([[int(marker.virtual_corners[0][0]),int(marker.virtual_corners[0][1])],
                                [int(marker.virtual_corners[1][0]),int(marker.virtual_corners[1][1])],
                                [int(marker.virtual_corners[2][0]),int(marker.virtual_corners[2][1])],
                                [int(marker.virtual_corners[3][0]),int(marker.virtual_corners[3][1])]
                                ]).reshape(1, -1, 2)
            frame = cv2.polylines(frame, points, isClosed=True, color=(0,0,256), thickness=2)
        else:
            points = np.array([[int(marker.corners[0][0]),int(marker.corners[0][1])],
                                [int(marker.corners[1][0]),int(marker.corners[1][1])],
                                [int(marker.corners[2][0]),int(marker.corners[2][1])],
                                [int(marker.corners[3][0]),int(marker.corners[3][1])]
                                ]).reshape(1, -1, 2)
            frame = cv2.polylines(frame, points, isClosed=True, color=(0,0,256), thickness=2)
    return frame

def load_checkpoint(checkpoint_file, model, optimizer, lr):
    print("=> Loading checkpoint")
    checkpoint = torch.load(checkpoint_file, map_location=config.DEVICE)
    model.load_state_dict(checkpoint["state_dict"])
    optimizer.load_state_dict(checkpoint["optimizer"])


def make_mesh(width, height,markers, lc, filename):
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

    gmsh.write(filename)
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
        if marker.id == config.ID_MAGNET:
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


def return_input_data(markers):
    # 学習データ用の白い背景画像を用意
    InputDataImage = np.zeros((config.HEIGHT, config.WIDTH, 3),dtype=np.uint8 ) + 255
    b1, g1, r1 = 255, 0, 0   # 青
    b2, g2, r2 = 0, 0, 255   # 赤
    b, g, r    = 0, 255, 0   # 緑
    for marker in markers:
        # print(marker.id)
        if marker.id == config.ID_MAGNET:
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
    return InputDataImage

