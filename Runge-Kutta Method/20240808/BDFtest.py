import numpy as np
import matplotlib.pyplot as plt

# Function definitions
def f(x):
    return np.exp(-x)

def g(x):
    return -np.exp(-x)

# Newton's method to solve u1 - dx*f(u1) - u0 = 0
def solve(u0, dx):
    y0 = u0
    err = 1e10
    while err > 1e-14:
        y1 = y0 - (y0 - u0 - dx * f(y0)) / (1 - dx * g(y0))
        err = abs(y1 - y0) / abs(y1)
        y0 = y1
    return y0

# Main program using Backward Euler method
def backward_euler():
    T = np.log(2)
    print("----------------------------------------")
    print(" log_2 dx         u        log_2(|error|) ")
    print("----------------------------------------")
    
    dx_vals = []
    u_vals = []
    error_vals = []
    
    for i in range(3, 13):
        dx = 2 ** (-i)
        u0 = 0.0
        x0 = 0.0
        
        while x0 + dx <= 1.0:
            u1 = solve(u0, dx)
            u1 = u0 + dx * f(u1)
            x1 = x0 + dx
            
            x0 = x1
            u0 = u1
        
        dx_vals.append(dx)
        u_vals.append(u0)
        error_vals.append(np.log2(abs(u0 - T)))
        
        print(f" {i:3d} {dx:15.7e} {u0:10.3e} {np.log2(abs(u0 - T)):10.3f}")
    
    print("----------------------------------------")
    
    # 結果のプロット
    plt.figure()
    plt.plot(np.log2(dx_vals), u_vals, label='u')
    plt.xlabel('log_2(dx)')
    plt.ylabel('u')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(np.log2(dx_vals), error_vals, label='log_2(|error|)')
    plt.xlabel('log_2(dx)')
    plt.ylabel('log_2(|error|)')
    plt.legend()
    plt.show()

# Run the main program
backward_euler()
