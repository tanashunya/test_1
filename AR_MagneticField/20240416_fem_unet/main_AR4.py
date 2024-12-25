import numpy as np
import cv2
from cv2 import aruco
# import math
# import gmsh
import Carray
import Cfunctions
from paraview.simple import *
import paraviewscript5 as para
import os
import time
# from datetime import datetime
import config
# import MyFunctions2
import MyFunctions3 as myfunc3
from MyFunctions3 import Marker
import torch
import torch.nn as nn
import torch.optim as optim
# from PIL import Image
# from torch.utils.data import Dataset, DataLoader
from torchvision.utils import save_image
# from tqdm import tqdm
import albumentations as A
from albumentations.pytorch import ToTensorV2
# import shutil
import matplotlib.pyplot as plt
import pandas as pd
# import csv
# import statistics
# import random
# import torchsummary
from torchvision.transforms import ToTensor, Resize
from MyGeneratorClass import *


# # 自動データ作成モード
# if config.MAKE_DATASET_AUTO:
#     myfunc3.make_dataset_auto(max_iteration=100000, fin_data_index=10000)
#     sys.exit()

class TimeMeasurement:
    def __init__(self):
        self.start_time = time.perf_counter()

    def print_time(self):
        self.time = time.perf_counter() - self.start_time
        print(self.time)
        

# DLモデル使う場合は、あらかじめモデルをロードしておく
if config.USE_DEEPLEARNING_MODEL:
    gen = Unet2(in_channels=3, out_channels=2, features=64).to(config.DEVICE)
    opt_gen = optim.Adam(gen.parameters(), lr=2e-4, betas=(0.5, 0.999))
    myfunc3.load_checkpoint(config.TAR_PATH_GEN, gen, opt_gen, 2e-4)

# 縦幅と横幅（解像度）
imgsize_x = config.WIDTH
imgsize_y = config.HEIGHT
height = float(imgsize_y)
width = float(imgsize_x)

# ウィンドウの設定
window_name = 'AR'
x_position = 1000
y_position = 0
# ウィンドウを作成
cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
cv2.resizeWindow(window_name, 1280, 960)
cv2.moveWindow(window_name, x_position, y_position)

# ウィンドウの設定
# window_name_input = 'input image'
# x_position = 0
# y_position = 0
# # ウィンドウを作成
# cv2.namedWindow(window_name_input, cv2.WINDOW_NORMAL)
# cv2.resizeWindow(window_name_input, 1280, 960)
# cv2.moveWindow(window_name_input, x_position, y_position)


# カメラの設定
cap = cv2.VideoCapture(0)  # 0はデフォルトのウェブカメラ
origin_y = 112
origin_x = 192

# ARマーカの設定
aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_50)
aruco_params = aruco.DetectorParameters_create()


# mesh size
lc = config.MESH_SIZE

# 現在のディレクトリを取得
current_dir = os.getcwd()

while True:
    loopstarttime = time.perf_counter()
    looptime = TimeMeasurement()
    captime = TimeMeasurement()
    markers = []
    ret, frame = cap.read()  # フレームのキャプチャ
    if config.WIDTH!=640 or config.HEIGHT!=480:
        frame = frame[origin_x:origin_x+config.WIDTH, origin_y:origin_y+config.HEIGHT]
    frame = cv2.imread("frametest.png") ###############################################################
    if config.ROTATE_180:
        frame = cv2.flip(frame, -1) # 180度回す
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)  # グレースケールへの変換
    corners, ids, rejected = aruco.detectMarkers(gray, aruco_dict, parameters=aruco_params)  # ARマーカの検出

    # cv2.imwrite("frametest.png", frame)

    # マーカが検出された場合、インスタンスを生成
    if ids is not None:
        num_mag = 0
        num_iron = 0
        for i in range(len(ids)):
            marker_id = ids[i][0]
            marker_corners = corners[i][0]
            m = Marker(marker_id, marker_corners, num_mag)
            markers.append(m)
            if myfunc3.is_sticking_out(m.virtual_corners):
                markers.pop(-1)
                print("はみ出しているオブジェクトを削除しました。")
            elif marker_id == config.ID_MAGNET:
                num_mag += 1
            elif marker_id == config.ID_IRON:
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
            
            if config.USE_DEEPLEARNING_MODEL:
                # Use Deep Learning model---------------------------------------------
                folder = "eval"
                input_image = myfunc3.return_input_data(markers)
                # cv2.imshow(window_name_input, input_image)
                # captime.print_time()
                dl_time = TimeMeasurement()
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
                    save_image(output, folder + f"/output_test.png")
                gen.train()
                # overlay_image = cv2.imread(folder + f"/output_test.png")
                # frame = cv2.imread(folder + f"/output_test.png")
                # frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)
                
                # dl_time.print_time()
                paraviewtimeDL = TimeMeasurement()

                # print(time.perf_counter() - start_time)
                if config.USE_PARAVIEW:
                    if config.DEVICE == "cuda":
                        output_numpy = output.squeeze().permute(1, 2, 0).cpu().numpy()
                    else:
                        output_numpy = output.squeeze().permute(1, 2, 0).numpy()
                    path_of_VTK = "DL_vector.vtk"
                    thinning_interval = int(config.WIDTH / config.NUMBER_OF_DIVISION)
                    # Bx = (output_numpy[1::thinning_interval, 1::thinning_interval, 0].flatten()) * 2.0*config.MAGNETIZATION - config.MAGNETIZATION
                    # By = (output_numpy[1::thinning_interval, 1::thinning_interval, 1].flatten()) * 2.0*config.MAGNETIZATION - config.MAGNETIZATION
                    Bx = (output_numpy[1::thinning_interval, 1::thinning_interval, 0].flatten())
                    By = (output_numpy[1::thinning_interval, 1::thinning_interval, 1].flatten())
                    
                    Bx, By = Bx[::-1], -1.0*By[::-1]
                    # print(Bx.shape) 
                    mycfunc = Cfunctions.Cfunctions()
                    mycfunc.write_vtk(
                        config.NUMBER_OF_DIVISION, config.NUMBER_OF_DIVISION,
                        Bx,By,
                        path_of_VTK
                    )
                    para.makevecpng_DL(current_dir, path_of_VTK)
                    # paraviewtimeDL.print_time()
                    overlay_image = cv2.imread("paraview_img_DL.png")
                    frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)


            else:
                # FEM ----------------------------------------------------------------
                meshtime = TimeMeasurement()

                meshname = "variablemesh2.msh1" 
                myfunc3.make_mesh(width, height, markers, lc, meshname)

                # meshtime.print_time()######
                femtime = TimeMeasurement()

                theta_of_magnets = np.array([])
                for marker in markers:
                    if marker.id == config.ID_MAGNET:
                        theta_of_magnets = np.append(theta_of_magnets, marker.thetarad)
                # print(theta_of_magnets)
                if(num_iron ==0 or config.NONLINER==False):
                    nonliner = False
                    # print(num_iron)
                else:
                    nonliner = True
                
                theta = 0.0
                my_carray = Carray.Carray(theta, 
                                        len(theta_of_magnets), 
                                        theta_of_magnets, 
                                        int(width), int(height), 
                                        config.MAKE_DATASET, i, 
                                        path_of_output_data,
                                        meshname)
                my_carray.read_mesh()
                my_carray.read_akimacsv()
                my_carray.analysis(nonliner)
                my_carray.write_vtk("vector")
                
                # femtime.print_time()######
                
                paraviewtime_fem = TimeMeasurement()
                # ベクトル or コンター描画--------------------------------------------------------
                if config.DRAW_VECTOR:
                    # frame = Myfunctions.draw_vectors(frame, imgsize_y) # Legacy ver.
                    para.makevecpng(current_dir)
                    
                    # paraviewtime_fem.print_time()######

                    overlay_image = cv2.imread("paraview_img.png")
                    frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)
                elif config.DRAW_CONTOUR:
                    para.makeconpng(current_dir)
                    overlay_image = cv2.imread("contour.png")
                    frame = cv2.addWeighted(frame, 1.0, overlay_image, 0.7, 0)


    cv2.imshow(window_name, frame)
    # looptime.print_time()######

    key = cv2.waitKey(1)
    if key == ord('q'):  # qを押すと終了
        if config.FRAME_SAVE_MODE:
            cv2.imwrite("frame.png", frame)
            print("frame.png saved successfully!")
        break
    
    # ループ終端


cap.release()  # カメラの解放
cv2.destroyAllWindows()  # ウィンドウを閉じる
