from scipy import fftpack
import numpy as np
import matplotlib.pyplot as plt

# 質量・減衰・剛性の集中定数を設定する
m1 = 0.5
m2 = 0.3
c1 = 0
c2 = 0
c12 = 0
k1 = 1.0e4
k2 = 1.0e2
k_12 = 45

Z = 1000000  # 回路インピーダンス
Amp = 1   # 外力の振幅
freq = 5  # 外力の周波数

# 運動方程式を関数として定義
def f(x, v):
    

    M = np.array([[m1, 0],              # 質量マトリクス
                 [0, m2]])
    C = np.array([[c1 + c12, -c12],       # 減衰マトリクス
                 [-c12, c2 + c12]])
    K = np.array([[k1 + k_12, -k_12],       # 剛性マトリクス
                 [-k_12, k2 + k_12]])
    M_inv = np.linalg.inv(M)            # 質量マトリクスの逆行列

    # 外力ベクトル
    F = np.array([[Amp * np.cos(2 * np.pi * freq * t)], 
                  [Amp * np.cos(2 * np.pi * freq * t)]])

    # 減衰強制振動モデル
    y = np.dot(M_inv, F) - np.dot((np.dot(C, M_inv)), v) - np.dot((np.dot(K, M_inv)), x)
    return y

# 解析条件
t0 = 0.0                            # 開始時間
t1 = 5.0                            # 終了時間
dt = 1e-4                          # 時間刻み

# 初期条件
x0 = np.array([[0.0], [0.0]])      # 初期変位[m]
v0 = np.array([[0.0], [0.0]])       # 初期速度[m/s]
x, v = x0, v0                       # 初期値に設定

t_axis = np.arange(t0, t1, dt)      # 時間軸
x1_sol = []                          # 初期化x(変位）配列
x2_sol = []                          # 初期化x(変位）配列
v1_sol = []                          # 初期化v(速度)配列
v2_sol = []                          # 初期化v(速度)配列

# 4次のRunge-Kutta法による数値解析
iteration = 0
for t in t_axis:
    x1_sol.append(x[0, 0])           # 1つ目の自由度の変位を結果として抽出
    x2_sol.append(x[1, 0])           # 2つ目の自由度の変位を結果として抽出
    v1_sol.append(v[0, 0])           # 1つ目の自由度の速度を結果として抽出
    v2_sol.append(v[1, 0])           # 2つ目の自由度の速度を結果として抽出
    print('iteration=', iteration,'time=', t)
    k11 = f(x, v) * dt
    k12 = v * dt
    k21 = f(x + k11 / 2, v + k12 / 2) * dt
    k22 = (v + k12 / 2) * dt
    k31 = f(x + k21 / 2, v + k22 / 2) * dt
    k32 = (v + k22 / 2) * dt
    k41 = f(x + k31, v + k32) * dt
    k42 = (v + k32) * dt
    v += (k11 + 2 * k21 + 2 * k31 + k41) / 6
    x += (k12 + 2 * k22 + 2 * k32 + k42) / 6
    iteration += 1

# 検証用の周波数分析
spec = fftpack.fft(x1_sol)                        # 時間波形をフーリエ変換してスペクトルにする
abs_spec = np.abs(spec / (len(spec)/2))          # 振幅を計算
frequency = np.linspace(0, 1/dt, len(x1_sol))     # 周波数軸を作成

# ここからグラフ描画
# フォントの種類とサイズを設定する。
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'Times New Roman'

# 目盛を内側にする。
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# グラフの上下左右に目盛線を付ける。
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax2 = fig.add_subplot(212)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')

# 軸のラベルを設定する。
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Displacement [m]')
ax2.set_xlabel('Frequency [Hz]')
ax2.set_ylabel('Displacement [m]')

ax2.set_xlim(0, 10)
ax2.set_yscale('log')

# データプロット
ax1.plot(t_axis, x1_sol, label='x1')
ax1.plot(t_axis, x2_sol, label='x2')
ax1.legend()  # ラベルを表示
ax2.plot(frequency, abs_spec)

fig.tight_layout()

# グラフを表示する。
plt.show()
plt.close()