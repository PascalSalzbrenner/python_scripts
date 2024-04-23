# script to plot the MAEs and RMSEs from EDDP hyperparameter tests as a function of cutoff radius and number of polynomials
# written by Pascal Salzbrenner, pts28@cam.ac.uk

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# read in the data
column_types = {'r': float, 'P': int, 'RMSE': float, 'MAE': float}
data = pd.read_csv('hyperparameter_test_average.dat', sep=' ', dtype=column_types)

# plot the MAEs
heatmap_data = data.pivot_table(index='P', columns='r', values='MAE')
heatmap = sns.heatmap(heatmap_data, cmap="magma_r")

plt.xlabel(r'r$_\text{c}$ [$\mathrm{\AA}$]', fontsize=12)
plt.ylabel('Number of Polynomials, M', fontsize=12)

colorbar = heatmap.collections[0].colorbar
colorbar.set_label('MAE [meV/atom]', fontsize=12)

x_ticks = np.arange(0, 551, step=50)
x_tick_labels = ['{:.1f}'.format(data['r'].unique()[tick]) for tick in x_ticks]

plt.xticks(ticks=x_ticks, labels=x_tick_labels, rotation=0)
plt.yticks(rotation=0)

plt.gca().invert_yaxis()

plt.subplots_adjust(top=0.92, bottom=0.12)

plt.savefig("mae_heatmap.png", dpi=300)

plt.clf()

# plot the RMSEs
heatmap_data = data.pivot_table(index='P', columns='r', values='RMSE')
heatmap = sns.heatmap(heatmap_data, cmap="magma_r")

plt.xlabel(r'r$_\text{c}$ [$\mathrm{\AA}$]', fontsize=12)
plt.ylabel('Number of Polynomials, M', fontsize=12)

colorbar = heatmap.collections[0].colorbar
colorbar.set_label('RMSE [meV/atom]', fontsize=12)

x_ticks = np.arange(0, 551, step=50)
x_tick_labels = ['{:.1f}'.format(data['r'].unique()[tick]) for tick in x_ticks]

plt.xticks(ticks=x_ticks, labels=x_tick_labels, rotation=0)
plt.yticks(rotation=0)

plt.gca().invert_yaxis()

plt.subplots_adjust(top=0.92, bottom=0.12)

plt.savefig("rmse_heatmap.png", dpi=300)

