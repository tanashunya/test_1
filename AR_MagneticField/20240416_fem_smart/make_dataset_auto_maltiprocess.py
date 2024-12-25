import config
import numpy as np
import cv2
import torch
import gmsh
import math
import random
import sys
import os
import Carray
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures
import multiprocessing
import functools
from MyFunctions3 import *
import time

def make_dataset_auto(max_iteration, initial_data_index, fin_data_index, mesh_filename, lock):
    
    mesh_size   = 1.5372
    num_mag     = 1
    num_iron    = 1
    base_size   = 30
    fluctuation = 3 # base_sizeに加える乱数値の範囲を決める値
    size        = base_size + random.uniform(-fluctuation, fluctuation) 
    
    for iter in range(max_iteration):
        # print(f"iter: {iter}")
        
        objects = []
        num_mag_actual = 0
        num_iron_actual = 0
        for i in range (num_mag):
            rand_obj = Random_Object(size=size)
            m = Marker(id=config.ID_MAGNET, corners=rand_obj.coordinates, num_mag=i)
            # m.print_info()
            objects.append(m)
            if is_sticking_out(m.virtual_corners):
                objects.pop(-1)
                # print("はみ出しているオブジェクトを削除しました。")
            else:
                num_mag_actual += 1
        
        for i in range (num_iron):
            rand_obj = Random_Object(size=size)
            m = Marker(id=2, corners=rand_obj.coordinates, num_mag=0)
            # m.print_info()
            objects.append(m)
            if is_sticking_out(m.virtual_corners):
                objects.pop(-1)
                # print("はみ出しているオブジェクトを削除しました。")
            else:
                num_iron_actual += 1
        
        if num_mag_actual==0:
            # ("磁石がないよ")
            pass
        elif not is_possible_to_make_mesh(objects):
            # print("メッシュを生成できないよ")
            pass
        elif num_mag_actual==num_mag and num_iron_actual==num_iron:
            i = initial_data_index
            with lock:
                while(os.path.exists(config.DIR_OF_INPUT_DATA+f"{i}.png")):
                    i += 1
                if i>fin_data_index:
                    exit()
                print(f"data_index = {i}")
                path_of_input_data = config.DIR_OF_INPUT_DATA+f"{i}.png"
                path_of_output_data = config.DIR_OF_OUTPUT_DATA+f"{i}.csv"
                make_input_data(objects, path_of_input_data)
            
            theta_of_magnets = np.array([])
            for object in objects:
                if object.id == config.ID_MAGNET:
                    # print(object.thetarad)
                    theta_of_magnets = np.append(theta_of_magnets, object.thetarad)
            make_mesh(float(config.WIDTH), float(config.HEIGHT), objects, mesh_size, filename=mesh_filename)            
            my_carray = Carray.Carray(0.0, 
                                    len(theta_of_magnets), 
                                    theta_of_magnets, 
                                    config.WIDTH, config.HEIGHT, 
                                    True, i, 
                                    path_of_output_data,
                                    mesh_filename)
            my_carray.fem()

            

if __name__=="__main__":
    start_time = time.time()
    max_workers = os.cpu_count()  # Get the number of CPU cores
    print(max_workers)
    mesh_filename1 = "test.msh1"
    # make_dataset_auto(max_iteration=100, initial_data_index=100001, fin_data_index=10005, mesh_filename=mesh_filename1)
    with multiprocessing.Manager() as manager:
        lock = manager.Lock()
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            processes = []
            for process_index in range(max_workers):
                p = executor.submit(
                    make_dataset_auto,
                    max_iteration=300000, initial_data_index=10000, fin_data_index=30000, mesh_filename=f'meshfiles/process{process_index}.msh1', lock=lock)
                processes.append(p)
    
    end_time = time.time()
    elapsed_time = end_time - start_time

    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)

    print(f"Total time taken: {int(hours)} hours {int(minutes)} minutes {int(seconds)} seconds")