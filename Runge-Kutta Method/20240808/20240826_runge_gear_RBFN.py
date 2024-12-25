import numpy as np
import matplotlib.pyplot as plt
from RBFN_model import Fm_RBFN, Phi_RBFN
import torch
import torch.nn as nn

# 定数の定義
m1 = 0.5  # kg
m2 = 0.3  # kg
k1 = 10000  # N/m
k2 = 1000  # N/m
k12 = 45  # N/m
c1 = 0  # ダンパー定数
c2 = 0  # ダンパー定数
c12 = 0  # ダンパー定数
i_f = 0.5
Z = 1000000  # 回路インピーダンス
Amp = 0.001   # 外力の振幅
freq = 1  # 外力の周波数


time_span = 1  # シミュレーション時間
dt = 0.0001  # 時間ステップ


global fm_model
global phi_model
fm_model = Fm_RBFN('F_m_final_model.pth')
phi_model = Phi_RBFN('Phi_sum_final_model.pth')


# 電磁力 Fm と鎖交磁束 Φ の関数定義
def Fm(z1, z2, i_f, i_s):
    return fm_model.predict([[z1- z2, i_f, i_s]])

def Phi(z1, z2, i_f, i_s):
    return phi_model.predict([[z1- z2, i_f, i_s]])

# 外力 Fext の関数定義
def Fext(t):
    return Amp * np.cos(2 * np.pi * freq * t)  # 正弦波

# 式(1), (2)の定義
def f_1_2(z1, z2, z1_dot, z2_dot, i_f, i_s, t):
    g1 = (z1_dot*(-c1 - c12) + z2_dot*(c12) + z1*(-k1 - k12) + z2*(k12) - Fm(z1, z2, i_f, i_s) - m1 * Fext(t)) / m1
    g2 = (z2_dot*(-c2 - c12) + z1_dot*(c12) + z2*(-k2 - k12) + z1*(k12) - Fm(z1, z2, i_f, i_s) - m2 * Fext(t)) / m2
    # g1 = (z1_dot*(-c1 - c12) + z2_dot*(c12) + z1*(-k1 - k12) + z2*(k12) - m1 * Fext(t)) / m1
    # g2 = (z2_dot*(-c2 - c12) + z1_dot*(c12) + z2*(-k2 - k12) + z1*(k12) - m2 * Fext(t)) / m2
    return g1, g2

# 式(3)の定義
def f_3(z1, z2, i_f, i_s):
    return -Z * i_s

# 前進オイラー法？
def forward_euler_step(phi_current, phi_prev, z1, z2, i_f, i_s, dt):
    phi_next = Phi(z1, z2, i_f, i_s)
    i_s = -(phi_next - phi_current)/(Z*dt)
    return phi_next, i_s

# ルンゲクッタ法による式(1)の数値解法
def runge_kutta_step(z1, z2, z1_dot, z2_dot, i_f, i_s, dt, t):
    k1_1 = z1_dot
    k1_2 = z2_dot
    l1_1, l1_2 = f_1_2(z1, 
                       z2, 
                       z1_dot, 
                       z2_dot, 
                       i_f, 
                       i_s, 
                       t)

    k2_1 = z1_dot + 0.5 * l1_1 * dt
    k2_2 = z2_dot + 0.5 * l1_2 * dt
    l2_1, l2_2 = f_1_2(z1 + 0.5 * k1_1 * dt, 
                       z2 + 0.5 * k1_2 * dt, 
                       z1_dot + 0.5 * l1_1 * dt, 
                       z2_dot + 0.5 * l1_2 * dt, 
                       i_f, 
                       i_s, 
                       t + 0.5 * dt)

    k3_1 = z1_dot + 0.5 * l2_1 * dt
    k3_2 = z2_dot + 0.5 * l2_2 * dt
    l3_1, l3_2 = f_1_2(z1 + 0.5 * k2_1 * dt, 
                       z2 + 0.5 * k2_2 * dt, 
                       z1_dot + 0.5 * l2_1 * dt, 
                       z2_dot + 0.5 * l2_2 * dt, 
                       i_f, 
                       i_s, 
                       t + 0.5 * dt)

    k4_1 = z1_dot + l3_1 * dt
    k4_2 = z2_dot + l3_2 * dt
    l4_1, l4_2 = f_1_2(z1 + k3_1 * dt, 
                       z2 + k3_2 * dt, 
                       z1_dot + l3_1 * dt, 
                       z2_dot + l3_2 * dt, 
                       i_f, 
                       i_s, 
                       t + dt)

    z1_next = z1 + (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) * dt / 6
    z2_next = z2 + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) * dt / 6
    z1_dot_next = z1_dot + (l1_1 + 2 * l2_1 + 2 * l3_1 + l4_1) * dt / 6
    z2_dot_next = z2_dot + (l1_2 + 2 * l2_2 + 2 * l3_2 + l4_2) * dt / 6

    return z1_next, z2_next, z1_dot_next, z2_dot_next

##################################################
if __name__ == "__main__":
    # 結果を格納する配列
    amplitudes_z1 = []
    amplitudes_z2 = []
    amplitudes_phi = []

    # 初期条件
    z1 = 0
    z2 = 0
    z1_dot = 0
    z2_dot = 0
    i_s = 0.00
    phi_current = Phi(z1, z2, i_f, i_s)
    phi_prev = phi_current  # 初期条件としてphi_prevをphi_currentと同じに設定

    # 時間領域
    t = np.arange(0, time_span, dt)
    num_steps = len(t)

    # 結果を格納する配列
    z1_vals = np.zeros(num_steps)
    z2_vals = np.zeros(num_steps)
    phi_vals = np.zeros(num_steps)
    i_s_vals = np.zeros(num_steps)

    # シミュレーションループ
    for i in range(num_steps):
        print(f"{i} / {num_steps}")
        z1, z2, z1_dot, z2_dot = runge_kutta_step(z1, z2, z1_dot, z2_dot, i_f, i_s, dt, t[i])
        phi_next, tmp = forward_euler_step(phi_current, phi_prev, z1, z2, i_f, i_s, dt)
        # 更新
        phi_prev = phi_current
        phi_current = phi_next
        
        z1_vals[i] = z1
        z2_vals[i] = z2
        phi_vals[i] = phi_current
        i_s_vals[i] = i_s
    # 結果のプロット
    plt.figure()
    plt.plot(t, z1_vals, label='z1')
    plt.plot(t, z2_vals, label='z2')
    plt.xlabel('Time [s]')
    plt.ylabel('Displacement [m]')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(t, phi_vals, label='Phi')
    plt.xlabel('Time [s]')
    plt.ylabel('Phi')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(t, i_s_vals, label='i_s')
    plt.xlabel('Time [s]')
    plt.ylabel('i_s')
    plt.legend()
    plt.show()