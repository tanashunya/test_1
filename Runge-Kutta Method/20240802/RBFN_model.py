import torch
import torch.nn as nn
import numpy as np

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
        y_np = data[:, 5]  # x-direction Magnetic Force
        x = torch.tensor(x_np, dtype=torch.float32)
        y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

        self.model = RBF(centers=x, out_dim=1)
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict(self, input_data):
        input_tensor = torch.tensor(input_data, dtype=torch.float32)
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
        y_np = data[:, 3] + data[:, 4]  # Phi_1 + Phi_2
        x = torch.tensor(x_np, dtype=torch.float32)
        y = torch.tensor(y_np, dtype=torch.float32).view(-1, 1)

        self.model = RBF(centers=x, out_dim=1)
        self.model.load_state_dict(torch.load(model_path))
        self.model.eval()

    def predict(self, input_data):
        input_tensor = torch.tensor(input_data, dtype=torch.float32)
        input_tensor[0][0] = abs(input_tensor[0][0])
        with torch.no_grad():
            result = self.model(input_tensor).item()
        return result