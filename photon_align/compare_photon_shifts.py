import numpy as np
import matplotlib.pyplot as plt

Bars = {};

IDs_new = []
means_new = []
IDs_old = []
means_old = []
colors = []
#with open("../bin/ler_final_tdc_params.txt","r") as f, open("../../bandsoft_tools/include/LER_TDC_shifts.txt") as g:
with open("../bin/test_tdc_v2.txt","r") as f, open("../../bandsoft_tools/include/LER_TDC_shifts.txt") as g:
    for line in f:
        parse = line.strip().split()
        sector = int(parse[0])
        layer = int(parse[1])
        comp = int(parse[2])
        
        Bars[ layer*100 + sector* 10 + comp ] = [float(parse[5]) ]

        IDs_new.append( layer*100 + sector*10 + comp )
        means_new.append( float(parse[5]) )
    for line in g:
        parse = line.strip().split()
        sector = int(parse[0])
        layer = int(parse[1])
        comp = int(parse[2])

        IDs_old.append( layer*100 + sector*10 + comp )
        means_old.append( float(parse[5]) )

        Bars[ layer*100 + sector* 10 + comp ].append( float(parse[5]) )
        if np.abs( Bars[ layer*100 + sector*10 + comp ][0] - Bars[ layer*100 + sector*10 + comp ][1] ) > 1:
            print( sector,layer,comp)


means_old = np.asarray(means_old)
means_new = np.asarray(means_new)
plt.scatter(IDs_new,means_old-means_new)
plt.grid(True)
plt.ylim([-1,1])
plt.xlabel('Bar ID [Layer*100 + Sector*10 + Component]',fontsize=13)
plt.ylabel('Photon Peak (Old - New) [ns]',fontsize=13)
plt.savefig('photon_diff_zoomed.pdf',bbox_inches='tight')
plt.show()
