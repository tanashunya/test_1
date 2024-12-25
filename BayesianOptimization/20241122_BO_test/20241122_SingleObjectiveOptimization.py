import optuna
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Objective function
def objective(x: float) -> float:
    return np.cos(x) * x**-1 * np.sin(x)

# Difine the search space
def objective_optuna(trial: optuna.trial.Trial) -> float:
    x = trial.suggest_float("x", -10, 10)
    return objective(x)

# Visualization function
def plot_optimization_process(i, ax, study):
    ax.clear()
    x = np.linspace(-10, 10, 100)
    y = objective(x)
    ax.plot(x, y, label='Objective Function')
    for t in study.trials[:i]:
        if t.state == optuna.trial.TrialState.COMPLETE:
            x_trial = t.params['x']
            y_trial = objective(x_trial)
            ax.scatter(x_trial, y_trial, color='red', zorder=5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    ax.set_title(f"Iteration {i+1}")

def main_function():
    # Define the optimization study
    study = optuna.create_study(direction="maximize", 
                                sampler=optuna.samplers.TPESampler(),
                                pruner=optuna.pruners.MedianPruner())
    study.optimize(objective_optuna, n_trials=30)

    # Create gif
    fig, ax = plt.subplots()
    ani = FuncAnimation(fig, plot_optimization_process, frames=len(study.trials), fargs=(ax, study))
    ani.save('optimization_process.gif', writer='pillow')
    plt.show()

    # Print the best parameter
    print(f"Best parameters: {study.best_params}")
    print(f"Best value: {study.best_value}" )


if __name__ == '__main__':
    main_function()