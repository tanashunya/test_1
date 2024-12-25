import torch
import torch.nn as nn
import torch.optim as optim
import os
import numpy as np
from PIL import Image
from torch.utils.data import Dataset, DataLoader
import albumentations as A
from albumentations.pytorch import ToTensorV2
import pandas as pd
import config


# Generratorクラス中で使うBlockクラスの定義
class Block(nn.Module):
    def __init__(self, in_channels, out_channels, down=True, act="relu", use_dropout=False):
        super(Block, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, 4, 2, 1, bias=False, padding_mode="reflect")
            if down
            else nn.ConvTranspose2d(in_channels, out_channels, 4, 2, 1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU() if act == "relu" else nn.LeakyReLU(0.2),
        )

        self.use_dropout = use_dropout
        self.dropout = nn.Dropout(0.5)
        self.down = down

    def forward(self, x):
        x = self.conv(x)
        return self.dropout(x) if self.use_dropout else x


# Generatorクラスの定義
class Generator(nn.Module):
    def __init__(self, in_channels=3, out_channels=2, features=64):
        super().__init__()
        self.initial_down = nn.Sequential(
            nn.Conv2d(in_channels, features, 4, 2, 1, padding_mode="reflect"),
            nn.LeakyReLU(0.2),                                                                 # ↓サイズ
        )                                                                                      # 128
        self.down1 = Block(features,   features*2, down=True, act="leaky", use_dropout=False)  # 64
        self.down2 = Block(features*2, features*4, down=True, act="leaky", use_dropout=False)  # 32
        self.down3 = Block(features*4, features*8, down=True, act="leaky", use_dropout=False)  # 16
        self.down4 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 8
        self.down5 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 4
        self.down6 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 2
        self.bottleneck = nn.Sequential(
            nn.Conv2d(features*8, features*8, 4, 2, 1), nn.ReLU(),                             # 1x1
        )
        self.up1 = Block(features*8,   features*8, down=False, act="relu", use_dropout=True)
        self.up2 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=True)
        self.up3 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=True)
        self.up4 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=False)
        self.up5 = Block(features*8*2, features*4, down=False, act="relu", use_dropout=False)
        self.up6 = Block(features*4*2, features*2, down=False, act="relu", use_dropout=False)
        self.up7 = Block(features*2*2, features,   down=False, act="relu", use_dropout=False)
        self.final_up = nn.Sequential(
            nn.ConvTranspose2d(features*2, out_channels, kernel_size=4, stride=2, padding=1),
            nn.Tanh(),
        )

    def forward(self, x):
        d1 = self.initial_down(x)
        d2 = self.down1(d1)
        d3 = self.down2(d2)
        d4 = self.down3(d3)
        d5 = self.down4(d4)
        d6 = self.down5(d5)
        d7 = self.down6(d6)
        bottleneck = self.bottleneck(d7)
        up1 = self.up1(bottleneck)
        up2 = self.up2(torch.cat([up1, d7], 1))
        up3 = self.up3(torch.cat([up2, d6], 1))
        up4 = self.up4(torch.cat([up3, d5], 1))
        up5 = self.up5(torch.cat([up4, d4], 1))
        up6 = self.up6(torch.cat([up5, d3], 1))
        up7 = self.up7(torch.cat([up6, d2], 1))
        return self.final_up(torch.cat([up7, d1], 1))

# Unet1クラスの定義
class Unet1(nn.Module):
    def __init__(self, in_channels=3, out_channels=2, features=64):
        super().__init__()
        self.initial_down = nn.Sequential(
            nn.Conv2d(in_channels, features, 4, 2, 1, padding_mode="reflect"),
            nn.LeakyReLU(0.2),                                                                 # ↓サイズ
        )                                                                                      # 128
        self.down1 = Block(features,   features*2, down=True, act="leaky", use_dropout=False)  # 64
        self.down2 = Block(features*2, features*4, down=True, act="leaky", use_dropout=False)  # 32
        self.down3 = Block(features*4, features*8, down=True, act="leaky", use_dropout=False)  # 16
        self.down4 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 8
        self.down5 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 4
        self.down6 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 2
        self.bottleneck = nn.Sequential(
            nn.Conv2d(features*8, features*8, 4, 2, 1), nn.ReLU(),                             # 1x1
        )
        self.up1 = Block(features*8,   features*8, down=False, act="relu", use_dropout=True)
        self.up2 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=True)
        self.up3 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=True)
        self.up4 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=False)
        self.up5 = Block(features*8*2, features*4, down=False, act="relu", use_dropout=False)
        self.up6 = Block(features*4*2, features*2, down=False, act="relu", use_dropout=False)
        self.up7 = Block(features*2*2, features,   down=False, act="relu", use_dropout=False)
        self.final_up = nn.Sequential(
            nn.ConvTranspose2d(features*2, out_channels, kernel_size=4, stride=2, padding=1),
            # nn.Tanh(),
        )

    def forward(self, x):
        d1 = self.initial_down(x)
        d2 = self.down1(d1)
        d3 = self.down2(d2)
        d4 = self.down3(d3)
        d5 = self.down4(d4)
        d6 = self.down5(d5)
        d7 = self.down6(d6)
        bottleneck = self.bottleneck(d7)
        up1 = self.up1(bottleneck)
        up2 = self.up2(torch.cat([up1, d7], 1))
        up3 = self.up3(torch.cat([up2, d6], 1))
        up4 = self.up4(torch.cat([up3, d5], 1))
        up5 = self.up5(torch.cat([up4, d4], 1))
        up6 = self.up6(torch.cat([up5, d3], 1))
        up7 = self.up7(torch.cat([up6, d2], 1))
        return self.final_up(torch.cat([up7, d1], 1))

# Unet2クラスの定義
class Unet2(nn.Module):
    def __init__(self, in_channels=3, out_channels=2, features=64):
        super().__init__()
        self.initial_down = nn.Sequential(
            nn.Conv2d(in_channels, features, 4, 2, 1, padding_mode="reflect"),
            nn.LeakyReLU(0.2),                                                                 # ↓サイズ
        )                                                                                      # 128
        self.down1 = Block(features,   features*2, down=True, act="leaky", use_dropout=False)  # 64
        self.down2 = Block(features*2, features*4, down=True, act="leaky", use_dropout=False)  # 32
        self.down3 = Block(features*4, features*8, down=True, act="leaky", use_dropout=False)  # 16
        self.down4 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 8
        self.down5 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 4
        self.down6 = Block(features*8, features*8, down=True, act="leaky", use_dropout=False)  # 2
        self.bottleneck = nn.Sequential(
            nn.Conv2d(features*8, features*8, 4, 2, 1), nn.ReLU(),                             # 1x1
        )
        self.up1 = Block(features*8,   features*8, down=False, act="relu", use_dropout=False)
        self.up2 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=False)
        self.up3 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=False)
        self.up4 = Block(features*8*2, features*8, down=False, act="relu", use_dropout=False)
        self.up5 = Block(features*8*2, features*4, down=False, act="relu", use_dropout=False)
        self.up6 = Block(features*4*2, features*2, down=False, act="relu", use_dropout=False)
        self.up7 = Block(features*2*2, features,   down=False, act="relu", use_dropout=False)
        self.final_up = nn.Sequential(
            nn.ConvTranspose2d(features*2, out_channels, kernel_size=4, stride=2, padding=1),
            # nn.Tanh(),
        )

    def forward(self, x):
        d1 = self.initial_down(x)
        d2 = self.down1(d1)
        d3 = self.down2(d2)
        d4 = self.down3(d3)
        d5 = self.down4(d4)
        d6 = self.down5(d5)
        d7 = self.down6(d6)
        bottleneck = self.bottleneck(d7)
        up1 = self.up1(bottleneck)
        up2 = self.up2(torch.cat([up1, d7], 1))
        up3 = self.up3(torch.cat([up2, d6], 1))
        up4 = self.up4(torch.cat([up3, d5], 1))
        up5 = self.up5(torch.cat([up4, d4], 1))
        up6 = self.up6(torch.cat([up5, d3], 1))
        up7 = self.up7(torch.cat([up6, d2], 1))
        return self.final_up(torch.cat([up7, d1], 1))
#############################################################################################################
# データセットの定義
# class MyDataset(Dataset):
#     def __init__(self, root_dir):
#         self.root_dir = root_dir
#         self.input_dir = self.root_dir + "/input"
#         self.target_dir = self.root_dir + "/target"
#         self.list_files = os.listdir(self.root_dir)
#         self.list_input_files = os.listdir(self.input_dir)
#         self.list_target_files = os.listdir(self.target_dir)
#         self.transform_input = A.Compose([
#             A.Resize(width=256, height=256),
#             A.Normalize(mean=[0.5, 0.5, 0.5], std=[0.5, 0.5, 0.5], max_pixel_value=255.0,),
#             ToTensorV2()
#         ])
#         self.transform_target = A.Compose([
#             A.Resize(width=256, height=256),
#             ToTensorV2()
#         ])

#     def __len__(self):
#         return len(self.list_input_files)

#     def __getitem__(self, index):
#         input_file = self.list_input_files[index]
#         input_number = int(input_file.split('.')[0])
#         input_path = os.path.join(self.input_dir, input_file)
#         input_image = np.array(Image.open(input_path))
        
#         target_file = f"{input_number}.csv"
#         target_path = os.path.join(self.target_dir, target_file)
#         target_data = ((pd.read_csv(target_path, header=None)).dropna(axis=1)).values
#         target_data_channel1 = target_data[:, ::2]  # even -> Bx
#         target_data_channel2 = target_data[:, 1::2] # odd  -> By
#         # target_data_channel3 = np.zeros((256,256))
#         # target_data = np.stack((target_data_channel1, target_data_channel2, target_data_channel3), axis=2)
#         target_data = np.stack((target_data_channel1, target_data_channel2), axis=2)
        
#         # print(input_path, target_path, )
        
#         input_image = self.transform_input(image=input_image)["image"].to(torch.float32) 
#         target_data = self.transform_target(image=target_data)["image"].to(torch.float32)
#         target_data = (target_data+config.MAGNETIZATION)/(config.MAGNETIZATION*2)
#         # target_data = target_data/config.MAGNETIZATION

#         # print(input_image)
        
#         return input_image, target_data
    