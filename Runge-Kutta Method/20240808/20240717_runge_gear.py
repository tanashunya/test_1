
import numpy as np
import matplotlib.pyplot as plt

# 定数の定義
m1 = 0.5  # kg
m2 = 0.3  # kg
k1 = 10000  # N/m
k2 = 10000  # N/m
k12 = 45  # N/m
c1 = 10  # ダンパー定数
c2 = 10  # ダンパー定数
c12 = 5  # ダンパー定数
Z = 1e6  # 回路インピーダンス
time_span = 10  # シミュレーション時間
dt = 0.01  # 時間ステップ

# 電磁力 Fm と鎖交磁束 Φ の関数定義
def Fm(z1, z2, if_current, is_current):
    # ここに電磁力の計算を実装
    return 1

def Phi(z1, z2, if_current, is_current):
    # ここに磁束の計算を実装
    return 1

# フーリエ変換のための関数定義
def fourier_transform(signal, dt):
    n = len(signal)
    freq = np.fft.fftfreq(n, d=dt)
    fft_signal = np.fft.fft(signal)
    return freq, np.abs(fft_signal)



# ルンゲクッタ法による式(1)の数値解法
def runge_kutta_step(z1, z2, z1_dot, z2_dot, if_current, is_current, dt):
    k1_1 = z1_dot
    k1_2 = z2_dot
    l1_1 = (-c1 * z1_dot - c12 * (z1_dot - z2_dot) - k1 * z1 - k12 * (z1 - z2) + Fm(z1, z2, if_current, is_current)) / m1
    l1_2 = (-c2 * z2_dot - c12 * (z2_dot - z1_dot) - k2 * z2 - k12 * (z2 - z1) + Fm(z1, z2, if_current, is_current)) / m2

    k2_1 = z1_dot + 0.5 * l1_1 * dt
    k2_2 = z2_dot + 0.5 * l1_2 * dt
    l2_1 = (-c1 * (z1_dot + 0.5 * l1_1 * dt) - c12 * (z1_dot + 0.5 * l1_1 * dt - z2_dot - 0.5 * l1_2 * dt) - k1 * (z1 + 0.5 * k1_1 * dt) - k12 * (z1 + 0.5 * k1_1 * dt - z2 - 0.5 * k1_2 * dt) + Fm(z1 + 0.5 * k1_1 * dt, z2 + 0.5 * k1_2 * dt, if_current, is_current)) / m1
    l2_2 = (-c2 * (z2_dot + 0.5 * l1_2 * dt) - c12 * (z2_dot + 0.5 * l1_2 * dt - z1_dot - 0.5 * l1_1 * dt) - k2 * (z2 + 0.5 * k1_2 * dt) - k12 * (z2 + 0.5 * k1_2 * dt - z1 - 0.5 * k1_1 * dt) + Fm(z1 + 0.5 * k1_1 * dt, z2 + 0.5 * k1_2 * dt, if_current, is_current)) / m2

    k3_1 = z1_dot + 0.5 * l2_1 * dt
    k3_2 = z2_dot + 0.5 * l2_2 * dt
    l3_1 = (-c1 * (z1_dot + 0.5 * l2_1 * dt) - c12 * (z1_dot + 0.5 * l2_1 * dt - z2_dot - 0.5 * l2_2 * dt) - k1 * (z1 + 0.5 * k2_1 * dt) - k12 * (z1 + 0.5 * k2_1 * dt - z2 - 0.5 * k2_2 * dt) + Fm(z1 + 0.5 * k2_1 * dt, z2 + 0.5 * k2_2 * dt, if_current, is_current)) / m1
    l3_2 = (-c2 * (z2_dot + 0.5 * l2_2 * dt) - c12 * (z2_dot + 0.5 * l2_2 * dt - z1_dot - 0.5 * l2_1 * dt) - k2 * (z2 + 0.5 * k2_2 * dt) - k12 * (z2 + 0.5 * k2_2 * dt - z1 - 0.5 * k2_1 * dt) + Fm(z1 + 0.5 * k2_1 * dt, z2 + 0.5 * k2_2 * dt, if_current, is_current)) / m2

    k4_1 = z1_dot + l3_1 * dt
    k4_2 = z2_dot + l3_2 * dt
    l4_1 = (-c1 * (z1_dot + l3_1 * dt) - c12 * (z1_dot + l3_1 * dt - z2_dot - l3_2 * dt) - k1 * (z1 + k3_1 * dt) - k12 * (z1 + k3_1 * dt - z2 - k3_2 * dt) + Fm(z1 + k3_1 * dt, z2 + k3_2 * dt, if_current, is_current)) / m1
    l4_2 = (-c2 * (z2_dot + l3_2 * dt) - c12 * (z2_dot + l3_2 * dt - z1_dot - l3_1 * dt) - k2 * (z2 + k3_2 * dt) - k12 * (z2 + k3_2 * dt - z1 - k3_1 * dt) + Fm(z1 + k3_1 * dt, z2 + k3_2 * dt, if_current, is_current)) / m2

    z1_next = z1 + (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) * dt / 6
    z2_next = z2 + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) * dt / 6
    z1_dot_next = z1_dot + (l1_1 + 2 * l2_1 + 2 * l3_1 + l4_1) * dt / 6
    z2_dot_next = z2_dot + (l1_2 + 2 * l2_2 + 2 * l3_2 + l4_2) * dt / 6

    return z1_next, z2_next, z1_dot_next, z2_dot_next

# Gearの2次後退差分法による式(2)の数値解法
def gear_second_order_backward_diff_step(z1, z2, if_current, is_current, dt):
    phi = Phi(z1, z2, if_current, is_current)
    is_next = is_current - dt * phi / Z
    return is_next
##################################################

# 結果を格納する配列
amplitudes_z1 = []
amplitudes_z2 = []
amplitudes_is = []

# 初期条件
z1 = 0
z2 = 0
z1_dot = 0
z2_dot = 0
if_current = 0
is_current = 0

# 時間領域
t = np.arange(0, time_span, dt)
num_steps = len(t)

# 結果を格納する配列
z1_vals = np.zeros(num_steps)
z2_vals = np.zeros(num_steps)
is_vals = np.zeros(num_steps)

# シミュレーションループ
for i in range(num_steps):
    z1, z2, z1_dot, z2_dot = runge_kutta_step(z1, z2, z1_dot, z2_dot, if_current, is_current, dt)
    is_current = gear_second_order_backward_diff_step(z1, z2, if_current, is_current, dt)
    
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
