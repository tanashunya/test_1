import torch
import torch.nn as nn
import numpy as np
import csv
from scipy.interpolate import RegularGridInterpolator

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


class Fm_RBFN:
    def __init__(self, model_path):
        data = np.loadtxt('sum.csv', delimiter=',', skiprows=1)
        x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
        x = torch.tensor(x_np, dtype=torch.float32)

        self.model = RBF(centers=x, out_dim=1)
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict(self, input_data):
        input_tensor = torch.tensor(input_data, dtype=torch.float32)
        input_tensor[0][0] *= 1000  # xを1000倍する　mmで学習しちゃったから
        if input_tensor[0][0] > 0:
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return result
        else:
            input_tensor[0][0] = abs(input_tensor[0][0])
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return -1 * result
        
class Fm_RBFN_2:
    def __init__(self, model_path):
        # ガウス基底を配置する位置をさらに細かく
        x_range = np.arange(0, 20.5, 0.5)
        i_f_range = np.arange(-2, 2.5, 0.5)
        i_s_range = np.arange(-2, 2.5, 0.5)
        X, I_F, I_S = np.meshgrid(x_range, i_f_range, i_s_range)
        centers = np.column_stack((X.ravel(), I_F.ravel(), I_S.ravel()))
        centers_tensor = torch.tensor(centers, dtype=torch.float32)

        self.model = RBF(centers=centers_tensor, out_dim=1)
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict(self, input_data):
        input_tensor = torch.tensor(input_data, dtype=torch.float32)
        # print(input_tensor.shape) #torch.Size([1, 3])
        input_tensor[0][0] *= 1000  # xを1000倍する　mmで学習しちゃったから
        if input_tensor[0][0] > 0:
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return result
        else:
            input_tensor[0][0] = abs(input_tensor[0][0])
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return -1 * result
        
class Phi_RBFN:
    def __init__(self, model_path):
        data = np.loadtxt('sum.csv', delimiter=',', skiprows=1)
        x_np = data[:, [0, 1, 2]]  # x, i_f, i_s
        x = torch.tensor(x_np, dtype=torch.float32)

        self.model = RBF(centers=x, out_dim=1)
        self.model.load_state_dict(torch.load(model_path), strict=False)
        self.model.eval()

    def predict(self, input_data):
        input_tensor = torch.tensor(input_data, dtype=torch.float32)
        input_tensor[0][0] *= 1000  # xを1000倍する　mmで学習しちゃったから
        if input_tensor[0][0] > 0:
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return result
        else:
            input_tensor[0][0] = abs(input_tensor[0][0])
            with torch.no_grad():
                result = self.model(input_tensor).item()
            return  result

class Fm_Linear:
    def __init__(self, csv_path):
        x, i_f, i_s, force_x = [], [], [], []
        with open(csv_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # ヘッダーをスキップ
            for row in reader:
                x.append(float(row[0]))
                i_f.append(float(row[1]))
                i_s.append(float(row[2]))
                force_x.append(float(row[5]))
        x, i_f, i_s, force_x = np.array(x), np.array(i_f), np.array(i_s), np.array(force_x)
        
        x /= 1000

        self.F_m_interpolator = RegularGridInterpolator(points=(np.unique(i_s), np.unique(i_f), np.unique(x)), 
                                                values=force_x.reshape(len(np.unique(i_s)), len(np.unique(i_f)), len(np.unique(x))),
                                                method='linear')

    def interpolate(self, x, i_f, i_s):
        if x > 0.02:
            return self.F_m_interpolator((i_s, i_f, 0.02))
        elif x < -0.02:
            return -1 * self.F_m_interpolator((i_s, i_f, 0.02))
        elif x > 0:
            return self.F_m_interpolator((i_s, i_f, x)) 
        else:
            return -1 * self.F_m_interpolator((i_s, i_f, abs(x)))
        

