import cv2

cap = None
for device_id in range(10):
    cap = cv2.VideoCapture(device_id)
    if cap.isOpened():
        print(f"デバイスID {device_id} のカメラに接続しました")
        break
    else:
        print(f"デバイスID {device_id} のカメラに接続できません")

if cap and cap.isOpened():
    ret, frame = cap.read()
    if not ret:
        print("フレームを取得できません")
    else:
        # フレーム処理
        pass

if cap:
    cap.release()