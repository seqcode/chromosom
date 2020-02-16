import sys
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pylab as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

file = sys.argv[1]
outbase = sys.argv[2]
colors = sys.argv[3]

cdict = {'red':   [[0.0,  0.0, 0.0],
                   [1.0,  1.0, 1.0]],
         'green': [[0.0,  0.0, 0.0],
                   [1.0,  0.0, 0.0]],
         'blue':  [[0.0,  0.0, 0.0],
                   [1.0,  0.0, 0.0]]}
newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)
rgba = newcmp(np.linspace(0, 1, 256))

with open(file) as f:
    ncols = len(f.readline().split('\t'))
print str(ncols)+' columns\n'

data = np.loadtxt(file, skiprows=1, usecols=range(1,ncols))

fig = plt.figure(figsize=[8,6]);
if colors == 1:
    ax = sns.heatmap(data, linewidth=0, cmap='Reds', vmin=2.0)
else:
    ax = sns.heatmap(data, linewidth=0, center=0)
outname = outbase + ".all.png"
fig.savefig(outname, dpi=300)

chr1data = data[1:906, 1:906]
fig = plt.figure(figsize=[8,6]);
if colors == 1:
    ax = sns.heatmap(chr1data, linewidth=0, cmap='Reds', vmin=2.0)
else:
    ax = sns.heatmap(chr1data, linewidth=0, center=0)
#plt.show()
outname = outbase + ".chr1.png"
fig.savefig(outname, dpi=300)

