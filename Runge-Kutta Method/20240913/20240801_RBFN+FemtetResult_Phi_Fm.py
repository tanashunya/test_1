EPOCHS    = 10000

SEED      = 42
save_model= False
best_model_path = 'best_model.pth'
final_model_path = 'final_model.pth'
load_model= True
load_model_path = 'final_model.pth'  


import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image
import cv2
import csv

class RBF(nn.Module):
    def __init__(self, centers, out_dim):
        super(RBF, self).__init__()
        self.centers = centers
        
        # betaの初期値計算
        # d_max = np.max([np.linalg.norm(c1 - c2) for c1 in centers.detach().numpy() for c2 in centers.detach().numpy()])
        # beta = 1 / (2 * (d_max / np.sqrt(2 * centers.size(0)))**2)
        
        # self.beta = nn.Parameter(torch.ones(centers.size(0)) * beta)
        self.beta = nn.Parameter(torch.ones(centers.size(0)))
        # self.beta = nn.Parameter(torch.ones(centers.size(0)) * beta, requires_grad=False) #betaずっと固定，全部同じ値
        self.linear = nn.Linear(centers.size(0), out_dim)

    def kernel_function(self, x, c):
        return torch.exp(-self.beta.view(1, -1) * torch.sum((x - c) ** 2, dim=1, keepdim=True))

    def forward(self, x):
        phi = torch.exp(-self.beta.view(1, -1) * torch.sum((x.unsqueeze(1) - self.centers) ** 2, dim=2))
        return self.linear(phi)
    
    def print_info(self):
        # print(f'beta: {self.beta}')
        pass
        

def RBFN_Fm():  
    LR        = 0.001

    # シード値を固定
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    # CSVファイルからデータを読み込む
    data = np.loadtxt('results_femtet/sum.csv', delimiter=',', skiprows=1)
    
    # 入力データ (x, i_f, i_s) と出力データ (Magnetic Force) を抽出
    x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
    y_np = data[:, 5]  # x-direction Magnetic Force

    x = torch.tensor(x_np, dtype=torch.float32)
    y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

    model = RBF(centers=x, out_dim=1) # データの座標にガウス基底を配置
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LR)

    best_loss = float('inf')
    best_model_state = None

    if load_model:
        model.load_state_dict(torch.load(best_model_path+'_Fm', weights_only=True))
        print('model loaded')
    else:
        loss_file_path = 'results_MagneticForce/loss_record.csv'
        
        with open(loss_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Epoch', 'Loss'])
        
        for epoch in range(EPOCHS):
            model.train()
            optimizer.zero_grad()
            outputs = model(x)
            loss = criterion(outputs, y)
            loss.backward()
            optimizer.step()
            
            with open(loss_file_path, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([epoch+1, loss.item()])
            
            if (epoch+1) % 100 == 0:
                model.print_info()
                print(f'Epoch [{epoch+1}/{EPOCHS}], Loss: {loss.item():.10f}')

            if save_model:
                if loss.item() < best_loss:
                    best_loss = loss.item()
                    best_model_state = model.state_dict()

                if best_model_state is not None:
                    torch.save(best_model_state, best_model_path+'_Fm')
                    print(f'epoch: {epoch+1}, best model saved at {best_model_path+"_Fm"}')
        
        torch.save(model.state_dict(), 'F_m_'+final_model_path)
        print(f'final model saved: {'F_m_'+final_model_path}')
        
    model.eval()
    with torch.no_grad():
        result = model(x)

    images = []
    x_range = np.linspace(x_np[:, 0].min(), x_np[:, 0].max(), 100)
    i_f_range = np.linspace(x_np[:, 1].min(), x_np[:, 1].max(), 100)
    i_s_range = np.linspace(x_np[:, 2].min(), x_np[:, 2].max(), 100)
    
    X, I_F = np.meshgrid(x_range, i_f_range)
    
    y_min, y_max = y_np.min(), y_np.max()

    for i, i_s in enumerate(i_s_range):
        XI_F_IS = np.column_stack((X.ravel(), I_F.ravel(), np.full_like(X.ravel(), i_s)))
        XI_F_IS_tensor = torch.tensor(XI_F_IS, dtype=torch.float32)
        with torch.no_grad():
            Z = model(XI_F_IS_tensor).numpy().reshape(X.shape)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # # 実際のデータをプロット
        # scatter = ax.scatter(x_np[:, 0], x_np[:, 1], y_np, c=x_np[:, 2], cmap='viridis', label='True')
        
        # 予測面をプロット
        ax.plot_surface(X, I_F, Z, alpha=0.7)
        
        ax.set_xlabel('x1-x2')
        ax.set_ylabel('i_f')
        ax.set_zlabel('Magnetic Force')
        ax.set_zlim(y_min, y_max)  # yの範囲を固定
        # plt.colorbar(scatter, label='i_s')
        plt.title(f'Magnetic Force Interpolation (i_s={i_s:.2f})')
        # plt.legend()
        plt.savefig(f'results_MagneticForce/RBFN_magnetic_force_{i}.png')
        plt.close(fig)
        
        images.append(Image.open(f'results_MagneticForce/RBFN_magnetic_force_{i}.png'))
    images[0].save('RBFN_magnetic_force_animation.gif', save_all=True, append_images=images[1:], duration=50, loop=0)


def RBFN_Phi_1():
    LR        = 0.00005
    EPOCHS    = 30000
    # シード値を固定
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    # CSVファイルからデータを読み込む
    data = np.loadtxt('results_femtet/sum.csv', delimiter=',', skiprows=1)
    
    # 入力データ (x, i_f, i_s) と出力データ (Phi_1) を抽出
    x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
    y_np = data[:, 3]  # Phi_1

    x = torch.tensor(x_np, dtype=torch.float32)
    y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

    model = RBF(centers=x, out_dim=1) # データの座標にガウス基底を配置
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LR)

    best_loss = float('inf')
    best_model_state = None

    if load_model:
        model.load_state_dict(torch.load('Phi_1_'+best_model_path, weights_only=True))
        print('model loaded')
    else:
        loss_file_path = 'results_Phi_1/loss_record.csv'
        
        with open(loss_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Epoch', 'Loss'])

        for epoch in range(EPOCHS):
            model.train()
            optimizer.zero_grad()
            outputs = model(x)
            loss = criterion(outputs, y)
            loss.backward()
            optimizer.step()

            with open(loss_file_path, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([epoch+1, loss.item()])
            
            if (epoch+1) % 100 == 0:
                model.print_info()
                print(f'Epoch [{epoch+1}/{EPOCHS}], Loss: {loss.item():.10f}')

            if save_model:
                if loss.item() < best_loss:
                    best_loss = loss.item()
                    best_model_state = model.state_dict()

                if best_model_state is not None:
                    torch.save(best_model_state, 'Phi_1_'+best_model_path)
                    print(f'epoch: {epoch+1}, best model saved at {best_model_path+"_Phi_1"}')
        
        torch.save(model.state_dict(), 'Phi_1_'+final_model_path)
        print(f'final model saved: {'Phi_1_'+final_model_path}')
        
    model.eval()
    with torch.no_grad():
        result = model(x)

    images = []
    x_range = np.linspace(x_np[:, 0].min(), x_np[:, 0].max(), 100)
    i_f_range = np.linspace(x_np[:, 1].min(), x_np[:, 1].max(), 100)
    i_s_range = np.linspace(x_np[:, 2].min(), x_np[:, 2].max(), 100)
    
    X, I_F = np.meshgrid(x_range, i_f_range)
    
    y_min, y_max = y_np.min(), y_np.max()
    
    for i, i_s in enumerate(i_s_range):
        XI_F_IS = np.column_stack((X.ravel(), I_F.ravel(), np.full_like(X.ravel(), i_s)))
        XI_F_IS_tensor = torch.tensor(XI_F_IS, dtype=torch.float32)
        with torch.no_grad():
            Z = model(XI_F_IS_tensor).numpy().reshape(X.shape)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # # 実際のデータをプロット
        # scatter = ax.scatter(x_np[:, 0], x_np[:, 1], y_np, c=x_np[:, 2], cmap='viridis', label='True')
        
        # 予測面をプロット
        ax.plot_surface(X, I_F, Z, alpha=0.7)
        
        ax.set_xlabel('x1-x2')
        ax.set_ylabel('i_f')
        ax.set_zlabel('Phi_1')
        ax.set_zlim(y_min, y_max)  # yの範囲を固定
        # plt.colorbar(scatter, label='i_s')
        plt.title(f'Phi_1 Interpolation (i_s={i_s:.2f})')
        # plt.legend()
        plt.savefig(f'results_Phi_1/RBFN_Phi_1_{i}.png')
        plt.close(fig)
        
        images.append(Image.open(f'results_Phi_1/RBFN_Phi_1_{i}.png'))
    images[0].save('RBFN_Phi_1_animation.gif', save_all=True, append_images=images[1:], duration=50, loop=0)

def RBFN_Phi_2():
    LR        = 0.00005
    EPOCHS    = 30000
    # シード値を固定
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    # CSVファイルからデータを読み込む
    data = np.loadtxt('results_femtet/sum.csv', delimiter=',', skiprows=1)
    
    # 入力データ (x, i_f, i_s) と出力データ (Phi_2) を抽出
    x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
    y_np = data[:, 4]  # Phi_2

    x = torch.tensor(x_np, dtype=torch.float32)
    y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

    model = RBF(centers=x, out_dim=1) # データの座標にガウス基底を配置
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LR)

    best_loss = float('inf')
    best_model_state = None

    if load_model:
        model.load_state_dict(torch.load('Phi_2_'+best_model_path, weights_only=True))
        print('model loaded')
    else:
        loss_file_path = 'results_Phi_2/loss_record.csv'
        
        with open(loss_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Epoch', 'Loss'])

        for epoch in range(EPOCHS):
            model.train()
            optimizer.zero_grad()
            outputs = model(x)
            loss = criterion(outputs, y)
            loss.backward()
            optimizer.step()

            with open(loss_file_path, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([epoch+1, loss.item()])
            
            if (epoch+1) % 100 == 0:
                model.print_info()
                print(f'Epoch [{epoch+1}/{EPOCHS}], Loss: {loss.item():.10f}')

            if save_model:
                if loss.item() < best_loss:
                    best_loss = loss.item()
                    best_model_state = model.state_dict()

                if best_model_state is not None:
                    torch.save(best_model_state, 'Phi_2_'+best_model_path)
                    print(f'epoch: {epoch+1}, best model saved at {'Phi_2_'+best_model_path}')
        
        torch.save(model.state_dict(), 'Phi_2_'+final_model_path)
        print(f'final model saved: {'Phi_2_'+final_model_path}')
        
    model.eval()
    with torch.no_grad():
        result = model(x)

    images = []
    x_range = np.linspace(x_np[:, 0].min(), x_np[:, 0].max(), 100)
    i_f_range = np.linspace(x_np[:, 1].min(), x_np[:, 1].max(), 100)
    i_s_range = np.linspace(x_np[:, 2].min(), x_np[:, 2].max(), 100)
    
    X, I_F = np.meshgrid(x_range, i_f_range)
    
    y_min, y_max = y_np.min(), y_np.max()
    
    for i, i_s in enumerate(i_s_range):
        XI_F_IS = np.column_stack((X.ravel(), I_F.ravel(), np.full_like(X.ravel(), i_s)))
        XI_F_IS_tensor = torch.tensor(XI_F_IS, dtype=torch.float32)
        with torch.no_grad():
            Z = model(XI_F_IS_tensor).numpy().reshape(X.shape)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # # 実際のデータをプロット
        # scatter = ax.scatter(x_np[:, 0], x_np[:, 1], y_np, c=x_np[:, 2], cmap='viridis', label='True')
        
        # 予測面をプロット
        ax.plot_surface(X, I_F, Z, alpha=0.7)
        
        ax.set_xlabel('x1-x2')
        ax.set_ylabel('i_f')
        ax.set_zlabel('Phi_2')
        ax.set_zlim(y_min, y_max)  # yの範囲を固定
        # plt.colorbar(scatter, label='i_s')
        plt.title(f'Phi_1 Interpolation (i_s={i_s:.2f})')
        # plt.legend()
        plt.savefig(f'results_Phi_2/RBFN_Phi_2_{i}.png')
        plt.close(fig)
        
        images.append(Image.open(f'results_Phi_2/RBFN_Phi_2_{i}.png'))
    images[0].save('RBFN_Phi_2_animation.gif', save_all=True, append_images=images[1:], duration=50, loop=0)


def RBFN_Phi_sum():
    LR        = 0.00005
    EPOCHS    = 30000
    # シード値を固定
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    # CSVファイルからデータを読み込む
    data = np.loadtxt('results_femtet/sum.csv', delimiter=',', skiprows=1)
    
    # 入力データ (x, i_f, i_s) と出力データ (Phi_1) を抽出
    x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
    y_np = data[:, 3] + data[:, 4]  # Phi_1 + Phi_2

    x = torch.tensor(x_np, dtype=torch.float32)
    y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

    model = RBF(centers=x, out_dim=1) # データの座標にガウス基底を配置
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LR)

    best_loss = float('inf')
    best_model_state = None

    if load_model:
        model.load_state_dict(torch.load('Phi_sum_'+load_model_path, weights_only=True))
        print('model loaded')
    else:
        loss_file_path = 'results_Phi_sum/loss_record.csv'
        
        with open(loss_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Epoch', 'Loss'])

        for epoch in range(EPOCHS):
            model.train()
            optimizer.zero_grad()
            outputs = model(x)
            loss = criterion(outputs, y)
            loss.backward()
            optimizer.step()
            
            with open(loss_file_path, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([epoch+1, loss.item()])
            
            if (epoch+1) % 100 == 0:
                model.print_info()
                print(f'Epoch [{epoch+1}/{EPOCHS}], Loss: {loss.item():.15f}')

            if save_model:
                if loss.item() < best_loss:
                    best_loss = loss.item()
                    best_model_state = model.state_dict()

                if best_model_state is not None:
                    torch.save(best_model_state, 'Phi_sum_'+best_model_path)
                    print(f'epoch: {epoch+1}, best model saved at {best_model_path+"_Phi_sum"}')
        
        torch.save(model.state_dict(), 'Phi_sum_'+final_model_path)
        print(f'final model saved: {'Phi_sum_'+final_model_path}')
        
    model.eval()
    with torch.no_grad():
        result = model(x)

    images = []
    x_range = np.linspace(x_np[:, 0].min(), x_np[:, 0].max(), 100)
    i_f_range = np.linspace(x_np[:, 1].min(), x_np[:, 1].max(), 100)
    i_s_range = np.linspace(x_np[:, 2].min(), x_np[:, 2].max(), 100)
    
    X, I_F = np.meshgrid(x_range, i_f_range)
    
    y_min, y_max = y_np.min(), y_np.max()
    
    for i, i_s in enumerate(i_s_range):
        XI_F_IS = np.column_stack((X.ravel(), I_F.ravel(), np.full_like(X.ravel(), i_s)))
        XI_F_IS_tensor = torch.tensor(XI_F_IS, dtype=torch.float32)
        with torch.no_grad():
            Z = model(XI_F_IS_tensor).numpy().reshape(X.shape)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # # 実際のデータをプロット
        # scatter = ax.scatter(x_np[:, 0], x_np[:, 1], y_np, c=x_np[:, 2], cmap='viridis', label='True')
        
        # 予測面をプロット
        ax.plot_surface(X, I_F, Z, alpha=0.7)
        
        ax.set_xlabel('x1-x2')
        ax.set_ylabel('i_f')
        ax.set_zlabel('Phi')
        ax.set_zlim(y_min, y_max)  # yの範囲を固定
        # plt.colorbar(scatter, label='i_s')
        plt.title(f'Phi Interpolation (i_s={i_s:.2f})')
        # plt.legend()
        plt.savefig(f'results_Phi_sum/RBFN_Phi_sum_{i}.png')
        plt.close(fig)
        
        images.append(Image.open(f'results_Phi_sum/RBFN_Phi_sum_{i}.png'))
    images[0].save('RBFN_Phi_sum_animation.gif', save_all=True, append_images=images[1:], duration=50, loop=0)


if __name__ == "__main__":
    # RBFN_Fm()
    # RBFN_Phi_1()
    # RBFN_Phi_2()
    RBFN_Phi_sum()
