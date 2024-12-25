from PIL import Image, ImageDraw, ImageFont
import numpy as np
import cv2

# 画像の読み込み
image_path = 'paraviewpng/fem.png'  # 画像のパスを指定
img = Image.open(image_path)
width, height = img.size  # 画像の幅と高さを取得

# 透過情報を持つ画像をPILで作成
img_pil = Image.new('RGBA', (width, height), (255, 255, 255, 0))
draw = ImageDraw.Draw(img_pil)

# テキストの描画
text = "Hello, Transparent!"
font_path = 'path/to/your/font.ttf'  # 使用するフォントのパスを指定
font_size = 36
font = ImageFont.truetype(font_path, font_size)
text_color = (255, 255, 255, 255)  # テキストの色 (RGBA)
text_width, text_height = draw.textsize(text, font=font)  # テキストのサイズを取得
text_position = ((width - text_width) // 2, (height - text_height) // 2)  # 中央に配置
draw.text(text_position, text, fill=text_color, font=font)

# PIL画像をOpenCV画像に変換
img_cv = cv2.cvtColor(np.array(img_pil), cv2.COLOR_RGBA2BGRA)

# OpenCVで表示または保存
cv2.imshow('Image with Transparent Text', img_cv)
cv2.waitKey(0)
# cv2.imwrite('output_image.png', img_cv)  # 画像を保存する場合
cv2.destroyAllWindows()
