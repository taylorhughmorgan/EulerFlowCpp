"""
@date: 2025-03-24
@author: Hugh Morgan
@description: parse the outputs from EulerFlowCpp
"""
import json
import numpy as np
import matplotlib.pyplot as plt

sol_name = "sedov_results.json"

#%% read the data
with open(sol_name, 'r') as fin:
    solution = json.load(fin)

sol_data = {}
for key, finname in solution['fields'].items():
    sol_data[key] = np.genfromtxt(finname, delimiter=',', missing_values='-nan(ind)', filling_values=np.NaN)

#%% plotting
id_max   = sol_data['pressure(Pa)'].shape[0]
test_idx = np.arange(0, id_max, 20)
rGrid__m = solution['grids']['grid(m)']
times    = solution['grids']['times(s)']
nPlots = len(sol_data)

fig, ax = plt.subplots(nrows=nPlots)

for i, (key, data) in enumerate(sol_data.items()):
    for idx in test_idx:
        time__ms = times[idx] * 1000
        if 'pressure' in key or 'energy' in key:
            ax[i].semilogy(rGrid__m, data[idx], label=f"time={time__ms:.2f}ms")
        else:
            ax[i].plot(rGrid__m, data[idx], label=f"time={time__ms:.2f}ms")
    ax[i].set_ylabel(key)
    ax[i].grid(True)

ax[nPlots-1].legend(fontsize=5)
ax[nPlots-1].set_xlabel("grid(m)")
# %%
