"""
https://drive.google.com/drive/folders/1qrdRp9YN4YZhpwfGdiGvUmB9iD0IfzQy
"""
import numpy as np
import matplotlib.pyplot as plt

# 微分方程式
def differential_equation(x, y):
    return x * y

# Euler法
def euler_method(y0, x0, xn, h):
    num_steps = int((xn - x0)/h) + 1
    x_values = np.linspace(x0, xn, num_steps)
    y_values = np.zeros(num_steps)
    y_values[0] = y0

    for i in range(1, num_steps):
        y_values[i] = y_values[i-1] + h*differential_equation(y_values[i-1], x_values[i-1])
    
    return x_values, y_values

# 2次ルンゲクッタ法
def runge_kutta_2nd_order(y0, x0, xn, h):
    num_steps = int((xn - x0) / h) + 1
    x_values = np.linspace(x0, xn, num_steps)
    y_values = np.zeros(num_steps)
    y_values[0] = y0

    for i in range(1, num_steps):
        k1 = h * differential_equation(y_values[i - 1], x_values[i - 1])
        k2 = h * differential_equation(y_values[i - 1] + 0.5 * k1, x_values[i - 1] + 0.5 * h)
        y_values[i] = y_values[i - 1] + k2

    return x_values, y_values

# 4次ルンゲクッタ法
def runge_kutta_4th_order(y0, x0, xn, h):
    num_steps = int((xn - x0) / h) + 1
    x_values = np.linspace(x0, xn, num_steps)
    y_values = np.zeros(num_steps)
    y_values[0] = y0

    for i in range(1, num_steps):
        k1 = h * differential_equation(y_values[i - 1], x_values[i - 1])
        k2 = h * differential_equation(y_values[i - 1] + 0.5 * k1, x_values[i - 1] + 0.5 * h)
        k3 = h * differential_equation(y_values[i - 1] + 0.5 * k2, x_values[i - 1] + 0.5 * h)
        k4 = h * differential_equation(y_values[i - 1] + k3, x_values[i - 1] + h)
        y_values[i] = y_values[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return x_values, y_values

# 厳密解の計算
def exact_solution(x):
    return np.exp(x**2 / 2)

# 初期条件と解の範囲
y0 = 1.0
x0 = 0.0
xn = 2.0
h = 0.1

# オイラー法での解
x_euler, y_euler = euler_method(y0, x0, xn, h)

# 2次ルンゲクッタ法での解
x_rk2, y_rk2 = runge_kutta_2nd_order(y0, x0, xn, h)

# 4次ルンゲクッタ法での解
x_rk4, y_rk4 = runge_kutta_4th_order(y0, x0, xn, h)

# 厳密解の計算
y_exact = exact_solution(np.linspace(x0, xn, int((xn - x0) / h) + 1))

# データをCSVファイルに保存
data = np.column_stack((x_euler, y_euler, y_rk2, y_rk4, y_exact))
header = 'x, Euler Method, 2nd Order Runge-Kutta, 4th Order Runge-Kutta, Exact Solution'
np.savetxt('numerical_solutions_exact.csv', data, delimiter=',', header=header, comments='')

# 結果のプロット
plt.figure(figsize=(10, 6))
plt.plot(x_euler, y_euler, label='Euler Method')
plt.plot(x_rk2, y_rk2, label='2nd Order Runge-Kutta')
plt.plot(x_rk4, y_rk4, label='4th Order Runge-Kutta')
plt.plot(x_euler, y_exact, label='Exact Solution', linestyle='dashed')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.show()