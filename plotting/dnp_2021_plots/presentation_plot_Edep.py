import numpy as np
import matplotlib.pyplot as plt

def Read(parse,Edep,SpB,S,B):
    # parse = ['Ptbin:', '1', 'EdepCut:', '6', 'S+B:', '3701', 'S:', '1831.86', 'B:', '1869.14']
    edep = int(parse[3])
    spb = float(parse[5])
    s = float(parse[7])
    b = float(parse[9])
    Edep.append(edep)
    SpB.append(spb)
    S.append(s)
    B.append(b)

fig, (aT, aM, aL) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2,2, 1]},figsize=(12,12))
fig.subplots_adjust(hspace=0,wspace=0)

aM.tick_params(labelbottom=False,direction='in')
aL.tick_params(direction='in')


#plt.figure(1,figsize=(9,6))
plt.sca(aL)
plt.xlim([3,22])
plt.ylim([1,25])
plt.xticks([5,10,15,20],fontsize=20)
plt.yticks([5,10,15,20],fontsize=20)
plt.ylabel(r'$\delta S/S$',fontsize=27,labelpad=12)
plt.xlabel('Minimum Energy deposition [MeVee]',fontsize=27)
plt.axvline(x=10,linestyle='--',color='black',linewidth=2,zorder=-1)

#plt.figure(2,figsize=(9,6))
plt.sca(aM)
plt.xlim([3,22])
plt.ylim([0,5])
plt.yticks([1,2,3,4],fontsize=20)
plt.ylabel('S/B',fontsize=27,labelpad=12)
plt.axvline(x=10,linestyle='--',color='black',linewidth=2,zorder=-1)

plt.sca(aT)
plt.xlim([3,22])
plt.ylim([10,7500])
plt.yscale('log')
plt.yticks([100,1000],fontsize=20)
plt.ylabel(r'Counts',fontsize=27,labelpad=12)
plt.text(4.5,30,'Background',color='red',fontsize=22)
plt.text(13,3000,'Signal',color='green',fontsize=22)
plt.axvline(x=10,linestyle='--',color='black',linewidth=2,zorder=-1)

#files = ['Edep_10pt2.txt','Edep_10pt4.txt','Edep_10pt6.txt']
files = ['Edep_10pt4_al1.txt','Edep_10pt4_al2.txt','Edep_10pt4_al3.txt']
#files = ['Edep_10pt4.txt']
labels=[r'$1.3 < \alpha_S < 1.4$',
        r'$1.4 < \alpha_S < 1.5$',
        r'$1.5 < \alpha_S < 1.6$'   ]
Ebeam = ['10.2','10.4','10.6']
mark = ['o','s','*']
col = ['blue','red','green']

Edep      = [[] for j in range(len(files))]
SpB       = [[] for j in range(len(files))]
S         = [[] for j in range(len(files))]
B         = [[] for j in range(len(files))]
for i in range(len(files)):
    Edep    [i] = [[] for j in range(3) ]
    SpB     [i] = [[] for j in range(3) ]
    S       [i] = [[] for j in range(3) ]
    B       [i] = [[] for j in range(3) ]

    
    j = -1
    with open(files[i],"r") as f:
        for line in f:
            parse = line.strip().split()

            pTbin = int(parse[1])
            Read(parse,Edep[i][pTbin],SpB[i][pTbin],S[i][pTbin],B[i][pTbin])
        
        #Edep[i][1] = np.asarray(Edep[i][1]) + (i-1)/100
        SpB[i][1] = np.asarray( SpB[i][1] )
        B[i][1] = np.asarray( B[i][1] )
        S[i][1] = np.asarray( S[i][1] )

        plt.sca(aT)
        plt.scatter(Edep[i][1],B[i][1],marker=mark[i],color='red')
        plt.scatter(Edep[i][1],S[i][1],marker=mark[i],color='green')

        
        plt.sca(aM)
        plt.scatter(Edep[i][1],S[i][1]/B[i][1],marker=mark[i],color='blue',label=labels[i])
        #plt.scatter(Edep[i][1],B[i][1],marker=mark[i],color=col[i])

        plt.sca(aL)
        plt.scatter(Edep[i][1],100*np.sqrt(SpB[i][1] + B[i][1])/S[i][1],marker=mark[i],color='blue')
        #plt.errorbar(DataSimX[i][0],DataSim[i][0],yerr=DataSimErr[i][0],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)





plt.sca(aL)

plt.sca(aM)

plt.legend(numpoints=1,loc=2,fontsize=16)

plt.sca(aT)
ax=plt.gca()
textstr = '\n'.join((
    r'$E_\mathrm{beam} = 10.4$ GeV/c',
    r'$p_e > %1.0f$ GeV/c' % (3, ),
    r'$Q^2 > %1.0f$ GeV$^2$/c$^4$' % (2, ),
    r'$W > %1.0f$ GeV/c$^2$' % (2, ),
    r'$\cos\theta_{nq} < %1.1f$' % (-0.8, ),
    r'$p_n > %1.2f$ GeV/c' % (0.25, ), 
    r"$W' > %1.1f$ GeV/c$^2$" % (1.8, ),
    r'$p_T < %1.1f$ GeV/c' % (0.1, ),
    ))
props = dict(boxstyle='round', facecolor='white', alpha=1,zorder=99)
plt.text(0.72, 0.95, textstr, transform=ax.transAxes, fontsize=15,verticalalignment='top', bbox=props)
plt.savefig('ratio_combo_firstbin.pdf',bbox_inches='tight')



plt.savefig('SBratio.pdf',bbox_inches='tight')
plt.show()
