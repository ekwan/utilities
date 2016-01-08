import numpy as np
import pylab as plt
#from matplotlib import cm

# read text from file
lines = np.recfromcsv("ts.csv")
filenames=[ lines[i][0] for i in range(0,len(lines)) ]
energies =[ lines[i][1] for i in range(0,len(lines)) ]
gradients=[ lines[i][2] for i in range(0,len(lines)) ]

fig = plt.figure()
ax1 = fig.add_subplot(111)
col = ax1.scatter(gradients,energies)

def onpick(event):
    print event.ind

fig.canvas.mpl_connect('pick_event', onpick)
plt.show()
