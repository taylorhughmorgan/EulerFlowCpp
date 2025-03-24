"""
@date: 2025-03-24
@author: Hugh Morgan
@description: parse the outputs from EulerFlowCpp
"""
import numpy as np
import matplotlib.pyplot as plt

solution = {'pressure' : "press_solution.csv",
            'density'  : "rho_solution.csv",
            'velocity' : 'vel_solution.csv'
}

#%% read the data
sol_data = {}
for key, finname in solution.items():
    sol_data[key] = np.genfromtxt(finname, delimiter=',', missing_values='-nan(ind)', filling_values=np.NaN)

#%% plotting
id_max = sol_data['pressure'].shape[0]
test_idx = np.arange(0, id_max, 20)
fig, ax = plt.subplots(nrows=len(solution))

for i, (key, data) in enumerate(sol_data.items()):
    for idx in test_idx:
        if key == 'pressure':
            ax[i].semilogy(data[idx], label=f"{idx}")
        else:
            ax[i].plot(data[idx], label=f"{idx}")
    ax[i].set_ylabel(key)
    ax[i].grid(True)
# %%
