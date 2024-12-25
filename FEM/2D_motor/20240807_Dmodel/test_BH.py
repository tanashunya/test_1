import pandas as pd
import numpy as np

input_file = '50A470/nyu_B^2.csv'
output_file = '50A470/H_B.csv'

# 入力CSVファイルの読み込み
data = pd.read_csv(input_file, header=None, names=['B^2', 'nu'])

# 磁束密度Bを計算
data['B'] = np.sqrt(data['B^2'])

# 磁界強度Hを計算
data['H'] = data['B'] * data['nu']

# 結果を新しいCSVファイルに出力
data[['H', 'B']].to_csv(output_file, index=False, header=False)

print(f"磁束密度Bと磁界強度Hの関係が{output_file}に出力されました。")