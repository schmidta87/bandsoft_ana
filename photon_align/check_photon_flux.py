import numpy as np
import matplotlib.pyplot as plt

IDs = []
peaks = []
colors = []
with open("../include/TDC_pass1v0_initbar.txt","r") as f:
    for line in f:
        parse = line.strip().split()
        sector = int(parse[0])
        layer = int(parse[1])
        comp = int(parse[2]) 

        photon_peak = float(parse[4])

        if( sector == 3 or sector == 4):
            colors.append('red')
        else:
            colors.append('blue')
        peaks.append(photon_peak)
        IDs.append( layer*100 + sector*10 + comp )


plt.scatter(IDs,peaks,color=colors,marker='o')
plt.ylim([0,175])
plt.grid(True)
plt.savefig('photon_peakheight.pdf')
plt.show()
