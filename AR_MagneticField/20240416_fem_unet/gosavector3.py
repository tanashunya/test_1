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
from MyGeneratorClass import *

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
    current_dir = os.getcwd()
    DATA_DIR = "/home/syubuntu1/デスクトップ/workspace/data/20231218/"
    DATA_DIR = "/media/syunya-linux/KIOXIA/20231218/"
    DATA_INDEX = 22974


    # モデルをロードしておく
    gen = Unet2(in_channels=3, out_channels=2, features=64).to(config.DEVICE)
    opt_gen = optim.Adam(gen.parameters(), lr=2e-4, betas=(0.5, 0.999))
    myfunc3.load_checkpoint(config.TAR_PATH_GEN, gen, opt_gen, 2e-4)


    folder = "eval"
        # input_image = myfunc3.return_input_data(markers)
    input_image = cv2.imread(DATA_DIR + f"input/{DATA_INDEX}.png")
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
    # output_DL = (output_DL-0.5)*2*config.MAGNETIZATION
    output_DL_x = -1.0 * output_DL[:, :, 0].flatten()
    output_DL_y = output_DL[:, :, 1].flatten()

        # FEM ----------------------------------------------------------------
        # meshname = "variablemesh2.msh1" 
        # myfunc3.make_mesh(width, height, markers, 1.5372, meshname)

        # femstart = time.perf_counter()
        # theta_of_magnets = np.array([])
        # for marker in markers:
        #     if marker.id == config.ID_MAGNET:
        #         theta_of_magnets = np.append(theta_of_magnets, marker.thetarad)
        # print(theta_of_magnets)
        # theta = 0.0
    path_of_output_data = DATA_DIR + f"/target/{DATA_INDEX}.csv"
        # my_carray = Carray.Carray(theta, 
        #                         len(theta_of_magnets), 
        #                         theta_of_magnets, 
        #                         int(width), int(height), 
        #                         True, 1, 
        #                         path_of_output_data,
        #                         meshname)
        # my_carray.fem()
        
        
    # CSVファイルを読み込む
    data_fem = np.genfromtxt(path_of_output_data, delimiter=',')
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
        f"vtk/output_fem_index_{DATA_INDEX}.vtk"
    )
    mycfunc.write_vtk(
        256, 256,
        output_DL_x, output_DL_y,
        f"vtk/output_DL_index_{DATA_INDEX}.vtk"
    )
    mycfunc.write_vtk(
        256, 256,
        B_gosa_x, B_gosa_y,
        f"vtk/B_gosa_index_{DATA_INDEX}.vtk"
    )
    # para.makevecpng_DL(current_dir, "B_gosa.png")




