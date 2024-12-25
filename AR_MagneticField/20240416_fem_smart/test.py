import concurrent.futures
import multiprocessing
import time

def worker_func(shared_var, lock):
    with lock:
        shared_var.value += 1
        print(f"Process ID: {multiprocessing.current_process().pid}, Shared Value: {shared_var.value}")

if __name__ == "__main__":
    start_time = time.time()  # プログラム開始時刻を記録

    # Managerを使って共有変数を作成
    with multiprocessing.Manager() as manager:
        shared_int = manager.Value('i', 0)
        lock = manager.Lock()  # Lockを作成

        # 共有変数を使う複数のプロセスを生成
        with concurrent.futures.ProcessPoolExecutor() as executor:
            processes = []
            for _ in range(5):
                p = executor.submit(worker_func, shared_int, lock)
                processes.append(p)

            # 終了待ち
            for p in processes:
                p.result()

    end_time = time.time()  # プログラム終了時刻を記録
    elapsed_time = end_time - start_time  # 経過時間を計算

    # 経過時間を「時」「分」「秒」に変換
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Total time taken: {int(hours)} hours {int(minutes)} minutes {int(seconds)} seconds")
