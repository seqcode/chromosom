import sys
import numpy as np
import seaborn as sns
import pandas
import matplotlib
import matplotlib.pylab as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

file = sys.argv[1]
outbase = sys.argv[2]

with open(file) as f:
    labels = f.readline().split('\t');
    ncols = len(labels);
print str(ncols)+' columns\n'

data = pandas.read_csv(file, sep='\s+', index_col=0);

fig = sns.clustermap(data, method='average', metric='euclidean', cmap="vlag", yticklabels=True, xticklabels=True, figsize=(15,15))
fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xmajorticklabels(), fontsize = 5)
fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_ymajorticklabels(), fontsize = 5)

outname = outbase + ".png"
fig.savefig(outname, dpi=300)
plt.show()

