import numpy as np
import matplotlib.pyplot as plt

def f(y):
    return y**3 - 3*y + 2  # 例: 非線形関数 f(y) = y^3 - 3y + 2

def df(y):
    return 3*y**2 - 3  # f(y) の導関数

def backward_euler_step(y_n, dt):
    NR_iter_max = 10
    NR_conv     = 1e-10

    # ニュートン法の初期値として、前ステップの値を使用
    y_next = y_n
    
    # ニュートン法の反復
    for iteration in range(NR_iter_max):
        F = y_next - y_n - dt * f(y_next)
        dF = 1 - dt * df(y_next)
        
        y_next_new = y_next - F / dF
        
        # 収束判定
        if abs(y_next_new - y_next) < NR_conv:
            print(f"反復回数 = {iteration + 1}")
            y_next = y_next_new
            break
        
        y_next = y_next_new
    return y_next

def gear_second_order_step(y_n, y_nminus1, dt):
    NR_iter_max = 10
    NR_conv     = 1e-10

    # ニュートン法の初期値として、前ステップの値を使用
    y_next = y_n
    
    # ニュートン法の反復
    for iteration in range(NR_iter_max):
        F = y_next - (4/3)*y_n + (1/3)*y_nminus1 - (2/3)*dt*f(y_next)
        dF = 1 - (2/3)*dt*df(y_next)
        
        y_next_new = y_next - F / dF
        
        # 収束判定
        if abs(y_next_new - y_next) < NR_conv:
            print(f"反復回数 = {iteration + 1}")
            y_next = y_next_new
            break
        
        y_next = y_next_new
    
    return y_next

# 初期条件とパラメータ設定
y0 = 0.0  # 初期値
t0 = 0.0  # 開始時間
T = 2.0   # 終了時間
dt = 0.01  # 時間刻み

# 時間ステップを進める
times = np.arange(t0, T + dt, dt)
ys_backward_euler = np.zeros(len(times))
ys_gear = np.zeros(len(times))
ys_backward_euler[0] = y0
ys_gear[0] = y0

# 1ステップ目は後退オイラー法で計算
ys_gear[1] = backward_euler_step(ys_gear[0], dt)

for i in range(1, len(times)):
    ys_backward_euler[i] = backward_euler_step(ys_backward_euler[i-1], dt)
    if i > 1:
        ys_gear[i] = gear_second_order_step(ys_gear[i-1], ys_gear[i-2], dt)

# 結果のプロット
plt.figure()
plt.plot(times, ys_backward_euler, label='Backward Euler')
plt.plot(times, ys_gear, label="Gear's Second Order")
plt.xlabel('Time [s]')
plt.ylabel('y')
plt.legend()
plt.show()

