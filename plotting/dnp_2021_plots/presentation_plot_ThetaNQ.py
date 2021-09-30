import numpy as np
import matplotlib.pyplot as plt

def Read(parse,T,SpB,B):
    T.append(float(parse[0]))
    SpB.append(float(parse[1]))
    B.append(float(parse[2]))


plt.figure(1,figsize=(12,8))
plt.xlim([143.5,180])
plt.ylim([0,1750])
plt.xticks([150,160,170,180],fontsize=24)
plt.yticks([0,500,1000,1500],fontsize=24)
plt.ylabel('Counts',fontsize=26,labelpad=12)
plt.xlabel(r'$\theta_\mathrm{nq}$ [Deg.]',fontsize=28)

files = ['ThetaNQ_10pt2_Edep10.txt']

ToF       = [[] for j in range(len(files))]
SpB       = [[] for j in range(len(files))]
B         = [[] for j in range(len(files))]
for i in range(len(files)):

    with open(files[i],"r") as f:
        for line in f:
            parse = line.strip().split()

            Read(parse,ToF[i],SpB[i],B[i])
        

        plt.figure(1)

        #ToF[i].insert(0,ToF[i][0]-0.8)
        #SpB[i].insert(0,0)
        #B[i].insert(0,0)
        #SpB[i].insert(len(ToF[i]),0)
        #B[i].insert(len(ToF[i]),0)
        #ToF[i].insert(len(ToF[i]),ToF[i][len(ToF[i])-1]+0.8)
        #ToF[i] = np.asarray(ToF[i]) + 0.4
        plt.step(ToF[i],SpB[i],color='blue',linewidth=3)
        plt.step(ToF[i],B[i],color='red',linewidth=3)





#plt.legend(numpoints=1,loc=4,fontsize=16)

plt.text(165,1250,"Off-time Events",color='blue',fontsize=24)
plt.text(160,1520,"Mixed Background",color='red',fontsize=24)

ax=plt.gca()
textstr = '\n'.join((
    r'$E_\mathrm{beam} = 10.2$ GeV/c',
    r'$p_e > %1.0f$ GeV/c' % (3, ),
    r'$Q^2 > %1.0f$ GeV$^2$/c$^4$' % (2, ),
    r'$W > %1.0f$ GeV/c$^2$' % (2, ),
    r'$E_\mathrm{dep} > %1.0f$ MeVee' % (10, ),
    r'$\cos\theta_{nq} < %1.1f$' % (-0.8, ),
    ))
props = dict(boxstyle='round', facecolor='white', alpha=1,zorder=99)
plt.text(0.34, 0.42, textstr, transform=ax.transAxes, fontsize=22,verticalalignment='top', bbox=props)
plt.savefig('thetanq_edep10.pdf',bbox_inches='tight')



plt.show()
