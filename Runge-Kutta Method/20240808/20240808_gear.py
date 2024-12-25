import numpy as np
import matplotlib.pyplot as plt

# 微分方程式の定義 (例: dy/dt = -y)
def f(t, y):
    return -y

# 初期条件
y0 = 1.0
t0 = 0.0
t_end = 10.0
dt = 0.001

# 時間ステップの数
n_steps = int((t_end - t0) / dt)

# 結果を格納する配列
t_values = np.linspace(t0, t_end, n_steps + 1)
y_values = np.zeros(n_steps + 1)
y_values[0] = y0

# Gearの2次後退差分法の初期化
y_prev = y0
y_curr = y0 + dt * f(t0, y0)  # 初期のEulerステップ

# Gearの2次後退差分法の適用
for i in range(1, n_steps):
    t = t_values[i]
    y_next = (4/3) * y_curr - (1/3) * y_prev + (2/3) * dt * f(t, y_curr)
    y_prev = y_curr
    y_curr = y_next
    y_values[i + 1] = y_next

# 結果のプロット
plt.plot(t_values, y_values, label='Gear 2nd Order BDF')
plt.xlabel('Time')
plt.ylabel('y')
plt.legend()
plt.show()