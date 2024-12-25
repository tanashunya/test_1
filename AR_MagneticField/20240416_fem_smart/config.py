import torch
# Trueなら１　Falseなら0
ROTATE_180               = 0 # webcamが奥から生えてるとき用
DROW_OBJECT_CONTOUR_MODE = 1 # 赤い枠つける

WIDTH                    = 256
HEIGHT                   = 256
VIRTUAL                  = 1 # 材料を長方形にするかどうか
FRAME_SAVE_MODE          = 1
MESH_SIZE                = 1.5372
# MESH_SIZE                = 13
MESH_SIZE                = 30
COLOR_BAR                = 1 # カラーバー描画するかどうか
COLOR_BAR_RAINBOW        = 0
DRAW_VECTOR              = 1
DRAW_CONTOUR             = 0
NONLINER                 = 1

MAKE_DATASET             = 0 # データセットを作るかどうか
# MAKE_DATASET_AUTO        = 0 # 自動でのデータセット作成
INITIAL_INDEX            = 1 # 何番から始める
DIR_OF_INPUT_DATA        = "../data/input/"
DIR_OF_OUTPUT_DATA       = "../data/target/"

USE_DEEPLEARNING_MODEL   = 0
SEED                     = 42
# TAR_PATH_GEN             = "/home/syunya-linux/デスクトップ/MyWorkspace/MyDeepLeraning/Pix2Pix/20231112_Unet/tar/20231113_15:16/epo500gen.pth.tar"
# TAR_PATH_GEN             = "/home/syunya-linux/デスクトップ/MyWorkspace/MyDeepLeraning/Pix2Pix/20231112_Unet/tar/20231211_08:30/epo500gen.pth.tar"
TAR_PATH_GEN             = "/home/syunya-linux/デスクトップ/MyWorkspace/MyDeepLeraning/Pix2Pix/20240119_Unet/tar/20240305_22:13/epo100gen.pth.tar"

# TAR_PATH_GEN             = "epo500gen.pth.tar"
# TAR_PATH_GEN             = "/home/syubuntu1/デスクトップ/epo500gen.pth.tar"
# DEVICE                   = "cuda" if torch.cuda.is_available() else "cpu"   
DEVICE                   = "cpu"   
USE_PARAVIEW             = 1
NUMBER_OF_DIVISION       = 32

MAGNETIZATION            = 1.25
ID_MAGNET                = 0
ID_IRON                  = 1

# ↓使えません
DARK_CONTOUR             = 0
