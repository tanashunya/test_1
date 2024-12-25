import numpy as np
import cv2
from cv2 import aruco
import Carray
import Cfunctions
from paraview.simple import *
import paraviewscript5 as para
import os
import time
import config
import MyFunctions3 as myfunc3
from MyFunctions3 import Marker
import torch
import torch.nn as nn
import torch.optim as optim
from torchvision.utils import save_image
import albumentations as A
from albumentations.pytorch import ToTensorV2
import matplotlib.pyplot as plt
import pandas as pd
# import csv
# import statistics
# import random
# import torchsummary
from torchvision.transforms import ToTensor, Resize
from MyGeneratorClass import Generator

class AR:
    def __init__(self):
        pass
    def set_config(self, width, height, window_name, x_position, y_position):
        self.width = width
        self.height = height
        self.window_name = window_name
        self.x_position = x_position
        self.y_position = y_position
    def set_width_and_height(self, width, height):
        self.width = width
        self.height = height
    def set_windowname(self, window_name):
        self.window_name = window_name
    def set_window_position(self, x_position, y_position):
        self.x_position = x_position
        self.y_position = y_position
    def create_window(self):
        cv2.namedWindow(self.window_name, cv2.WINDOW_NORMAL)
    def move_window(self,  x_position, y_position):
        cv2.moveWindow(self.window_name, x_position, y_position)




if __name__ == "__main__":
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

    # ARインスタンス生成
    ar1 = AR()
    ar1.set_config(config.WIDTH,
                   config.HEIGHT,
                   'test',
                   300,
                   0)
    ar1.create_window()

    # ARマーカの設定
    aruco_dict = aruco.Dictionary_get(aruco.DICT_4X4_50)
    aruco_params = aruco.DetectorParameters_create()

    # カメラの設定
    cap = cv2.VideoCapture(0)  # 0はデフォルトのウェブカメラ
    origin_y = 112
    origin_x = 192

    # 現在のディレクトリを取得
    current_dir = os.getcwd()


    markers = []
    ret, frame = cap.read()  # フレームのキャプチャ
    if config.WIDTH!=640 or config.HEIGHT!=480:
        frame = frame[origin_x:origin_x+config.WIDTH, origin_y:origin_y+config.HEIGHT]
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
            print(marker_id)
            marker_corners = corners[i][0]
            # print(marker_corners)
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
    
    if num_mag==0:
        print("磁石がありません")
    elif myfunc3.is_possible_to_make_mesh(markers) is False and config.VIRTUAL:
        print("メッシュを生成できませんでした")
    else:
        if config.DROW_OBJECT_CONTOUR_MODE:
            frame = myfunc3.draw_obj_contour(frame, markers)

        # if config.MAKE_DATASET:
        #     i = config.INITIAL_INDEX
        #     while(os.path.exists(config.DIR_OF_INPUT_DATA+f"{i}.png")):
        #         i += 1
        #     print(f"data_index = {i}")
        #     path_of_input_data = config.DIR_OF_INPUT_DATA+f"{i}.png"
        #     path_of_output_data = config.DIR_OF_OUTPUT_DATA+f"{i}.csv"
        #     myfunc3.make_input_data(markers, path_of_input_data)
        path_of_input_data = "input_test.png"
    
        # if config.USE_DEEPLEARNING_MODEL:
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
            save_image(output_contatinated, folder + f"/output_test.png")
        gen.train()
        output_DL = (output.squeeze().permute(1, 2, 0).numpy()).astype(np.float32)
        # print(output_DL.shape)
        output_DL = (output_DL-0.5)*2*config.MAGNETIZATION
        output_DL_x = -1.0 * output_DL[:, :, 0].flatten()
        output_DL_y = output_DL[:, :, 1].flatten()

        # FEM ----------------------------------------------------------------
        meshname = "variablemesh2.msh1" 
        myfunc3.make_mesh(width, height, markers, 1.5372, meshname)

        femstart = time.perf_counter()
        theta_of_magnets = np.array([])
        for marker in markers:
            if marker.id == config.ID_MAGNET:
                theta_of_magnets = np.append(theta_of_magnets, marker.thetarad)
        # print(theta_of_magnets)
        theta = 0.0
        path_of_output_data = "csv/fem_256x256.csv"
        my_carray = Carray.Carray(theta, 
                                len(theta_of_magnets), 
                                theta_of_magnets, 
                                int(width), int(height), 
                                True, 1, 
                                path_of_output_data,
                                meshname)
        my_carray.fem()
        
        
        # CSVファイルを読み込む
        # data_fem = np.genfromtxt(path_of_output_data, delimiter=',')
        data_fem = ((pd.read_csv(path_of_output_data, header=None)).dropna(axis=1)).values
        data_channel1 = data_fem[:, ::2]  # even -> Bx
        # print(data_channel1.shape)
        data_channel2 = data_fem[:, 1::2] # odd  -> By
        output_fem = np.stack((data_channel1, data_channel2), axis=2)
        # data_channel1 = data_channel1[::-1]
        data_channel1 = data_channel1.reshape(256*256)
        # data_channel2 = data_channel2[::-1]
        data_channel2 = -1.0 * data_channel2.reshape(256*256)
        # データの形状を確認する
        # print(output_fem.shape) 
        ###############################################################


        # 3次元配列を2次元に変形する
        # output_DL_reshaped = output_DL.reshape((256*256, 2))
        # print(output_DL_reshaped.shape)
        # output_fem_reshaped = output_fem.reshape((256*256, 2))

        # print(f"Unet最大値:{np.max(output_DL_reshaped)}")
        # print(f"FEM最大値:{np.max(output_fem_reshaped)}")

        # hikizan
        # B_gosa = output_fem_reshaped - output_DL_reshaped
        # print(f"誤差最大値:{np.max(B_gosa)}")
        # print(f"誤差最小値:{np.min(B_gosa)}")

        # # B_gosa_reshaped = B_gosa.reshape((256 * 256, 2))

        # # データをCSVファイルとして保存する
        # np.savetxt('csv/output_DL.csv', output_DL_reshaped, delimiter=',')
        # np.savetxt('csv/output_fem.csv', output_fem_reshaped, delimiter=',')
        # np.savetxt('csv/B_gosa.csv', B_gosa, delimiter=',')

        # 
        thinning_interval = 1
        # B_gosa_x = (B_gosa[1::thinning_interval, 1::thinning_interval, 0].flatten())
        # B_gosa_y = (B_gosa[1::thinning_interval, 1::thinning_interval, 1].flatten())
        # B_gosa_x = B_gosa[0::thinning_interval, 0]
        B_gosa_x = data_channel1 - output_DL_x
        # B_gosa_y = B_gosa[0::thinning_interval, 1]
        B_gosa_y = data_channel2 - output_DL_y

        mycfunc = Cfunctions.Cfunctions()
        mycfunc.write_vtk(
            256, 256,
            data_channel1,
            data_channel2,
            "vtk/output_fem.vtk"
        )
        mycfunc.write_vtk(
            256, 256,
            output_DL_x, output_DL_y,
            "vtk/output_DL.vtk"
        )
        mycfunc.write_vtk(
            256, 256,
            B_gosa_x, B_gosa_y,
            "vtk/B_gosa.vtk"
        )
        # para.makevecpng_DL(current_dir, "B_gosa.png")




