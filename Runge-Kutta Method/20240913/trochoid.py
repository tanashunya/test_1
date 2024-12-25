import numpy as np
import matplotlib.pyplot as plt

# トロコイドの計算
def trochoid(a, b, theta):
    x = a * theta - b * np.sin(theta)
    y = a - b * np.cos(theta)
    return x, y

# トロコイドの微分の計算
def trochoid_derivative(a, b, theta):
    dx = a - b * np.cos(theta)
    dy = b * np.sin(theta)
    return dx, dy

# パラメータ設定
a = 1.0  # 波の高さ
b = 0.7  # 波の長さ
theta = np.linspace(0, 8 * np.pi, 1000)  # 角度の範囲

# トロコイドの計算
x, y = trochoid(a, b, theta)

# トロコイドの微分の計算
dx, dy = trochoid_derivative(a, b, theta)

# グラフ描画
plt.figure(figsize=(10, 6))

# トロコイドのグラフ
plt.subplot(2, 1, 1)
plt.plot(x, y, label='Trochoid')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Trochoid Wave Simulation')
plt.legend()
plt.grid(True)

# トロコイドの微分のグラフ
plt.subplot(2, 1, 2)
plt.plot(theta, dx, label='dx/dtheta')
plt.plot(theta, dy, label='dy/dtheta')
plt.xlabel('Theta')
plt.ylabel('Derivative')
plt.title('Derivative of Trochoid')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

