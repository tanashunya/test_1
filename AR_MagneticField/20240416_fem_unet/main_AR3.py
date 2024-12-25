import numpy as np
import cv2
from cv2 import aruco
import math
import gmsh
import Carray
import Cfunctions
from paraview.simple import *
import paraviewscript5 as para
import os
import time
from datetime import datetime
import config
# import MyFunctions2
import MyFunctions3 as myfunc3
import torch
import torch.nn as nn
import torch.optim as optim
from PIL import Image
from torch.utils.data import Dataset, DataLoader
from torchvision.utils import save_image
from tqdm import tqdm
import albumentations as A
from albumentations.pytorch import ToTensorV2
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import csv
import statistics
import random
import torchsummary
from torchvision.transforms import ToTensor, Resize
from MyGeneratorClass import Generator

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



# DLモデル使う場合は、あらかじめモデルをロードしておく
if config.USE_DEEPLEARNING_MODEL:
    gen = Generator(in_channels=3, out_channels=2, features=64).to(config.DEVICE)
    opt_gen = optim.Adam(gen.parameters(), lr=2e-4, betas=(0.5, 0.999))
    myfunc3.load_checkpoint(config.TAR_PATH_GEN, gen, opt_gen, 2e-4)

# 縦幅と横幅（解像度）
imgsize_x = config.WIDTH
imgsize_y = config.HEIGHT
height = float(imgsize_y)
width = float(imgsize_x)

# ウィンドウの設定
window_name = 'press q to exit'
x_position = 300
y_position = 0
# ウィンドウを作成
cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
# cv2.resizeWindow(window_name, 1280, 960)
cv2.moveWindow(window_name, x_position, y_position)

# カメラの設定
cap = cv2.VideoCapture(0)  # 0はデフォルトのウェブカメラ

# ARマーカの設定
aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_50)
aruco_params = aruco.DetectorParameters_create()


# mesh size
lc = config.MESH_SIZE

# 現在のディレクトリを取得
current_dir = os.getcwd()

counter=0
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
            elif marker_id == config.ID_MAGNET:
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
                myfunc3.make_input_data(markers, path_of_input_data)
            
            if config.USE_DEEPLEARNING_MODEL:
                # Use Deep Learning model---------------------------------------------
                folder = "eval"
                input_image = myfunc3.return_input_data(markers)
                transform_input = A.Compose([
                            A.Resize(width=256, height=256),
                            A.Normalize(mean=[0.5, 0.5, 0.5], std=[0.5, 0.5, 0.5], max_pixel_value=255.0,),
                            ToTensorV2()
                        ])
                input_image = transform_input(image=input_image)["image"].to(torch.float32) 
                input_image = input_image.unsqueeze(0)
                input_image = input_image.to(config.DEVICE)

                gen.eval()
                with torch.no_grad():
                    output = gen(input_image)
                    new_channel = 0.5 * torch.ones_like(output[:, :1, :, :])  # 0.5で初期化
                    output_contatinated = torch.cat([output, new_channel], dim=1) # 新しいチャネルを連結
                    # print(output.size())
                    # PyTorchテンソルをnumpy配列に変換し、次元の順序を変更する
                    # frame = cv2.cvtColor((output_contatinated.squeeze().permute(1, 2, 0).numpy() * 255).astype(np.uint8), cv2.COLOR_RGB2BGR)
                    # cv2.imwrite("eval/frame.png", frame)
                    save_image(input_image , folder + f"/input_test.png")
                    # save_image(output, folder + f"/output_test.png")
                gen.train()
                # overlay_image = cv2.imread(folder + f"/output_test.png")
                # frame = cv2.imread(folder + f"/output_test.png")
                # frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)
                
                # print(time.perf_counter() - start_time)
                if config.USE_PARAVIEW:
                    output_numpy = output.squeeze().permute(1, 2, 0).numpy()
                    path_of_VTK = "DL_vector.vtk"
                    thinning_interval = int(config.WIDTH / config.NUMBER_OF_DIVISION)
                    Bx = (output_numpy[1::thinning_interval, 1::thinning_interval, 0].flatten()) * 2.0*config.MAGNETIZATION - config.MAGNETIZATION
                    By = (output_numpy[1::thinning_interval, 1::thinning_interval, 1].flatten()) * 2.0*config.MAGNETIZATION - config.MAGNETIZATION
                    
                    Bx, By = Bx[::-1], -1.0*By[::-1]
                    # print(Bx.shape) 
                    mycfunc = Cfunctions.Cfunctions()
                    mycfunc.write_vtk(
                        config.NUMBER_OF_DIVISION, config.NUMBER_OF_DIVISION,
                        Bx,By,
                        path_of_VTK
                    )
                    para.makevecpng_DL(current_dir, path_of_VTK)
                    start_time = time.perf_counter()
                    overlay_image = cv2.imread("paraview_img_DL.png")
                    frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)


            else:
                # FEM ----------------------------------------------------------------
                myfunc3.make_mesh(width, height, markers, lc)
                femstart = time.perf_counter()
                
                theta_of_magnets = np.array([])
                for marker in markers:
                    if marker.id == config.ID_MAGNET:
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
    # print(time.perf_counter() - start_time)

    if cv2.waitKey(1) == ord('q'):  # qを押すと終了
        if config.FRAME_SAVE_MODE:
            cv2.imwrite("frame.png", frame)
            print("frame.png saved successfully!")
        break
    
    # ループ終端
    print(time.perf_counter() - loopstarttime)

    
    # counter+=1
    # if(counter>3):break


cap.release()  # カメラの解放
cv2.destroyAllWindows()  # ウィンドウを閉じる
