import numpy as np
import matplotlib.pyplot as plt
from RBFN_model import Fm_RBFN, Phi_RBFN
import torch
import torch.nn as nn

# 定数の定義
m1 = 0.5  # kg
m2 = 0.3  # kg
k1 = 10000  # N/m
k2 = 100  # N/m
k12 = 45  # N/m
c1 = 0  # ダンパー定数
c2 = 0  # ダンパー定数
c12 = 0  # ダンパー定数
time_span = 10  # シミュレーション時間
dt = 0.001  # 時間ステップ

Z = 1000000  # 回路インピーダンス
Amp = 1   # 外力の振幅
freq = 5  # 外力の周波数

global fm_model
global phi_model
fm_model = Fm_RBFN('F_m_final_model.pth')
phi_model = Phi_RBFN('Phi_sum_final_model.pth')


# 電磁力 Fm と鎖交磁束 Φ の関数定義
def Fm(z1, z2, if_current, is_current):
    return fm_model.predict([[z1- z2, if_current, is_current]])

def Phi(z1, z2, if_current, is_current):
    # return phi_model.predict([[z1- z2, if_current, is_current]])
    return 0.01

# フーリエ変換のための関数定義
def fourier_transform(signal, dt):
    n = len(signal)
    freq = np.fft.fftfreq(n, d=dt)
    fft_signal = np.fft.fft(signal)
    return freq, np.abs(fft_signal)

# 外力 Fext の関数定義
def Fext(t):
    return Amp * np.sin(2 * np.pi * freq * t)  # 正弦波

# 式(1)の定義
def differential_equation_1(z1, z2, z1_dot, z2_dot, if_current, is_current, t):
    g1 = (-c1 * z1_dot - c12 * (z1_dot - z2_dot) - k1 * z1 - k12 * (z1 - z2) + Fm(z1, z2, if_current, is_current) - m1 * Fext(t)) / m1
    g2 = (-c2 * z2_dot - c12 * (z2_dot - z1_dot) - k2 * z2 - k12 * (z2 - z1) + Fm(z1, z2, if_current, is_current) - m2 * Fext(t)) / m2
    return g1, g2

# 式(2)の定義
def differential_equation_2(z1, z2, if_current, is_current):
    phi = Phi(z1, z2, if_current, is_current)
    return (-phi / Z)
    

# ルンゲクッタ法による式(1)の数値解法
def runge_kutta_step(z1, z2, z1_dot, z2_dot, if_current, is_current, dt, t):
    k1_1 = z1_dot
    k1_2 = z2_dot
    l1_1, l1_2 = differential_equation_1(z1, 
                                         z2, 
                                         z1_dot, 
                                         z2_dot, 
                                         if_current, 
                                         is_current, 
                                         t)

    k2_1 = z1_dot + 0.5 * l1_1 * dt
    k2_2 = z2_dot + 0.5 * l1_2 * dt
    l2_1, l2_2 = differential_equation_1(z1 + 0.5 * k1_1 * dt, 
                                         z2 + 0.5 * k1_2 * dt, 
                                         z1_dot + 0.5 * l1_1 * dt, 
                                         z2_dot + 0.5 * l1_2 * dt, 
                                         if_current, 
                                         is_current, 
                                         t + 0.5 * dt)

    k3_1 = z1_dot + 0.5 * l2_1 * dt
    k3_2 = z2_dot + 0.5 * l2_2 * dt
    l3_1, l3_2 = differential_equation_1(z1 + 0.5 * k2_1 * dt, 
                                         z2 + 0.5 * k2_2 * dt, 
                                         z1_dot + 0.5 * l2_1 * dt, 
                                         z2_dot + 0.5 * l2_2 * dt, 
                                         if_current, 
                                         is_current, 
                                         t + 0.5 * dt)

    k4_1 = z1_dot + l3_1 * dt
    k4_2 = z2_dot + l3_2 * dt
    l4_1, l4_2 = differential_equation_1(z1 + k3_1 * dt, 
                                         z2 + k3_2 * dt, 
                                         z1_dot + l3_1 * dt, 
                                         z2_dot + l3_2 * dt, 
                                         if_current, 
                                         is_current, 
                                         t + dt)

    z1_next = z1 + (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) * dt / 6
    z2_next = z2 + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) * dt / 6
    z1_dot_next = z1_dot + (l1_1 + 2 * l2_1 + 2 * l3_1 + l4_1) * dt / 6
    z2_dot_next = z2_dot + (l1_2 + 2 * l2_2 + 2 * l3_2 + l4_2) * dt / 6

    return z1_next, z2_next, z1_dot_next, z2_dot_next

# Gearの2次後退差分法による式(2)の数値解法
def gear_second_order_backward_diff_step(z1, z2, if_current, is_current, is_prev, dt):
    is_next = (4/3) * is_current - (1/3) * is_prev + (2/3) * dt * differential_equation_2(z1, z2, if_current, is_current)
    return is_next
##################################################
if __name__ == "__main__":
    

    # 結果を格納する配列
    amplitudes_z1 = []
    amplitudes_z2 = []
    amplitudes_is = []

    # 初期条件
    z1 = 0
    z2 = 0
    z1_dot = 0
    z2_dot = 0
    if_current = 0.5
    is_current = 0
    is_prev = 0

    # 時間領域
    t = np.arange(0, time_span, dt)
    num_steps = len(t)

    # 結果を格納する配列
    z1_vals = np.zeros(num_steps)
    z2_vals = np.zeros(num_steps)
    is_vals = np.zeros(num_steps)

    # シミュレーションループ
    for i in range(num_steps):
        print(f"{i} / {num_steps}")
        z1, z2, z1_dot, z2_dot = runge_kutta_step(z1, z2, z1_dot, z2_dot, if_current, is_current, dt, t[i])
        is_next = gear_second_order_backward_diff_step(z1, z2, if_current, is_current, is_prev, dt)
        
        # 更新
        is_prev = is_current
        is_current = is_next
        
        z1_vals[i] = z1
        z2_vals[i] = z2
        is_vals[i] = is_current

    # 結果のプロット
    plt.figure()
    plt.plot(t, z1_vals, label='z1')
    plt.plot(t, z2_vals, label='z2')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(t, is_vals, label='is')
    plt.xlabel('Time [s]')
    plt.ylabel('Current [A]')
    plt.legend()
    plt.show()