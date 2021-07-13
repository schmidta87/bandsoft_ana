import matplotlib.pyplot as plt
import numpy as np


i = 0
Ratio = {}
with open("must_plot_bincontent.txt","r") as f:
    for line in f:
        parse = line.strip().split()
        
        alpha_s     = float(parse[0])
        this_ratio  = float(parse[1])

        if alpha_s in Ratio.keys():
            Ratio[alpha_s].append( this_ratio )
        else:
            Ratio[alpha_s] = [ this_ratio ]

Xps = [0.3,0.4,0.5,0.6,0.7]
for key in Ratio:
    if key == 1.25: continue
    plt.errorbar( Xps, Ratio[key] ,label = r'$\alpha_S$ = '+str(key),linestyle='None',marker='o')

plt.ylim([0.7,3])
plt.xlabel(r"x'",fontsize=20)
plt.ylabel(r"$R_{tag}/R_{inc} / (R_{x'=0.3})$",fontsize=20)
plt.grid(True)
plt.legend(numpoints=1,loc='best')
plt.savefig('ratio_tag_inc_xp03_inverted.pdf',bbox_inches='tight')
plt.show()
