import numpy as np
import scipy.io
import scipy as sp
import matplotlib as mpl
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from plotly.offline import plot
import pandas as pd
import math
import orca
import seaborn as sns

# Get path of the data to be used
inpath = '/Users/joanduprez/Desktop/W/Research/UR1-EA4712/LFP/Data_LFP emotions/' 
outpath = '/Users/joanduprez/Desktop/W/Research/UR1-EA4712/LFP/Data_LFP emotions/results/'

## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## BEHAVIOR
# load the data

LFPdat = pd.read_csv(inpath + '/behavior_data.txt', sep="\t", header=0)

ax = [None] * (2 + 2)
fig = plt.figure(figsize=(5, 3.5))

gs = mpl.gridspec.GridSpec(nrows=1,
                           ncols=1,
                           figure=fig,
                           wspace=0.25, hspace=0.25
                           )
ax[0] = fig.add_subplot(gs[0, 0])
sns.set_palette('Paired')
sns.violinplot(x="emotion", y="accuracy", hue='task', data=LFPdat, split=True)

sns.swarmplot(x="emotion", y="accuracy", hue='task', data=LFPdat, dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[0].set_ylabel('Accuracy rate')
#ax[0].legend().set_visible(False)
ax[0].set_xlabel('')

leg = plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01),
            borderaxespad=0, frameon=False, fontsize=8, ncol =2)

leg.get_texts()[0].set_text('Gender task')
leg.get_texts()[1].set_text('Emotion task')
leg.legendHandles[2].set_visible(False)
leg.legendHandles[3].set_visible(False)

leg.get_texts()[2].set_text('')
leg.get_texts()[3].set_text('')

ax[0].xaxis.set_ticklabels(['Neutral', 'Fear'])

fig.savefig(outpath + '/accuracy_violin_split.png',
            format='png', dpi=1000)


## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# Power
# load the data

LFPdat = pd.read_csv(inpath + '/statfile.txt', sep="\t", header=0)

# Plot

ax = [None] * (2 + 2)
fig = plt.figure(figsize=(10, 7))

gs = mpl.gridspec.GridSpec(nrows=2,
                           ncols=2,
                           figure=fig,
                           wspace=0.25, hspace=0.25
                           )
ax[0] = fig.add_subplot(gs[0, 0])
sns.set_palette('Paired')
sns.violinplot(x="emo", y="pwr", hue='task', data=LFPdat[LFPdat.freq == 'delta'], split=True)

sns.swarmplot(x="emo", y="pwr",hue='task', data=LFPdat[LFPdat.freq == 'delta'], dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[0].set_ylabel('delta power (dB)')
#ax[0].legend().set_visible(False)
ax[0].set_xlabel('')
ax[0].set_ylim([-12, 12])
# ax[0].collections[0].set_facecolor('tab:blue')
# ax[0].collections[1].set_facecolor('lightsteelblue')
# ax[0].collections[3].set_facecolor('tab:blue')
# ax[0].collections[4].set_facecolor('lightsteelblue')
leg = plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01),
            borderaxespad=0, frameon=False, fontsize=8, ncol =2)

leg.get_texts()[0].set_text('Gender task')
leg.get_texts()[1].set_text('Emotion task')
leg.legendHandles[2].set_visible(False)
leg.legendHandles[3].set_visible(False)

leg.get_texts()[2].set_text('')
leg.get_texts()[3].set_text('')

ax[0].xaxis.set_ticklabels(['Neutral', 'Fear'])
ax[0].text(-0.95, 14, 'A', fontsize=16, fontweight='bold')


ax[1] = fig.add_subplot(gs[0, 1])
sns.violinplot(x="emo", y="pwr", hue='task', data=LFPdat[LFPdat.freq == 'alpha'], split=True)
sns.swarmplot(x="emo", y="pwr",hue='task', data=LFPdat[LFPdat.freq == 'alpha'], dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[1].set_ylabel('alpha power (dB)')
ax[1].legend().set_visible(False)
ax[1].set_xlabel('')
ax[1].set_ylim([-12, 12])
# ax[1].collections[0].set_facecolor('tab:blue')
# ax[1].collections[1].set_facecolor('lightsteelblue')
# ax[1].collections[3].set_facecolor('tab:blue')
# ax[1].collections[4].set_facecolor('lightsteelblue')
ax[1].xaxis.set_ticklabels(['Neutral', 'Fear'])
ax[1].text(-0.95, 14, 'B', fontsize=16, fontweight='bold')

ax[2] = fig.add_subplot(gs[1, 0])
sns.violinplot(x="emo", y="pwr", hue='task', data=LFPdat[LFPdat.freq == 'beta'], split=True)
sns.swarmplot(x="emo", y="pwr",hue='task', data=LFPdat[LFPdat.freq == 'beta'], dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[2].set_ylabel('beta power (dB)')
ax[2].set_xlabel('')
ax[2].legend().set_visible(False)
ax[2].set_ylim([-12, 12])
# ax[2].collections[0].set_facecolor('tab:blue')
# ax[2].collections[1].set_facecolor('lightsteelblue')
# ax[2].collections[3].set_facecolor('tab:blue')
# ax[2].collections[4].set_facecolor('lightsteelblue')
ax[2].xaxis.set_ticklabels(['Neutral', 'Fear'])
ax[2].text(-0.95, 14, 'C', fontsize=16, fontweight='bold')

fig.savefig(outpath + '/power_violin_split.png',
            format='png', dpi=1000)


## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------
# Same with ITPC
# load the data

LFPdat = pd.read_csv(inpath + 'statfile_itpc.txt', sep="\t", header=0)

# Plot

ax = [None] * (2 + 2)
fig = plt.figure(figsize=(10, 3.5))

gs = mpl.gridspec.GridSpec(nrows=1,
                           ncols=2,
                           figure=fig,
                           wspace=0.25, hspace=0.25
                           )
ax[0] = fig.add_subplot(gs[0, 0])
sns.set_palette('Paired')
sns.violinplot(x="emo", y="pwr", hue='task', data=LFPdat[LFPdat.freq == 'delta'], split=True)

sns.swarmplot(x="emo", y="pwr",hue='task', data=LFPdat[LFPdat.freq == 'delta'], dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[0].set_ylabel('delta ITPCz')
#ax[0].legend().set_visible(False)
ax[0].set_xlabel('')
#ax[0].set_ylim([-12, 12])
# ax[0].collections[0].set_facecolor('tab:blue')
# ax[0].collections[1].set_facecolor('lightsteelblue')
# ax[0].collections[3].set_facecolor('tab:blue')
# ax[0].collections[4].set_facecolor('lightsteelblue')
leg = plt.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01),
            borderaxespad=0, frameon=False, fontsize=8, ncol =2)

leg.get_texts()[0].set_text('Gender task')
leg.get_texts()[1].set_text('Emotion task')
leg.legendHandles[2].set_visible(False)
leg.legendHandles[3].set_visible(False)

leg.get_texts()[2].set_text('')
leg.get_texts()[3].set_text('')

ax[0].xaxis.set_ticklabels(['Neutral', 'Fear'])
ax[0].set(ylim =(-200, 800))
ax[0].text(-0.95, 850, 'A', fontsize=16, fontweight='bold')


ax[1] = fig.add_subplot(gs[0, 1])
sns.violinplot(x="emo", y="pwr", hue='task', data=LFPdat[LFPdat.freq == 'theta'], split=True)
sns.swarmplot(x="emo", y="pwr",hue='task', data=LFPdat[LFPdat.freq == 'theta'], dodge=True,
              color="black", size=3, edgecolor='white', linewidth=0.5, alpha=0.5)
ax[1].set_ylabel('theta  ITPCz ')
ax[1].legend().set_visible(False)
ax[1].set_xlabel('')
#ax[1].set_ylim([-12, 12])
# ax[1].collections[0].set_facecolor('tab:blue')
# ax[1].collections[1].set_facecolor('lightsteelblue')
# ax[1].collections[3].set_facecolor('tab:blue')
# ax[1].collections[4].set_facecolor('lightsteelblue')
ax[1].xaxis.set_ticklabels(['Neutral', 'Fear'])
ax[1].set(ylim =(-200, 800))
ax[1].text(-0.75, 850, 'B', fontsize=16, fontweight='bold')

fig.savefig(outpath + '/itpc_violin_split.png',
            format='png', dpi=1000)

