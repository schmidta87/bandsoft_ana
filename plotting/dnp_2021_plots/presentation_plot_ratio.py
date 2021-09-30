import numpy as np
import matplotlib.pyplot as plt

def Read(line,X,Y,E):
    parse = line.strip().split()
    x = float(parse[0])
    y = float(parse[1])
    e = float(parse[2])
    X.append(x)
    Y.append(y)
    E.append(e)

def WeightedAverage( AvgX, Avg, AvgE, X, Y, E ):

    for pT_bin in range(len(X[0])):
        for as_bin in range(len(X[0][pT_bin])):
            for xp_bin in range(len(X[0][pT_bin][as_bin])):
                avg = 0
                wei = 0
                cnt = 0

                # sum over the three beam energies:
                for en_bin in range(len(X)):

                    if E[en_bin][pT_bin][as_bin][xp_bin] == 0: continue
                    avg += Y[en_bin][pT_bin][as_bin][xp_bin] / pow(E[en_bin][pT_bin][as_bin][xp_bin],2)
                    wei += 1./pow(E[en_bin][pT_bin][as_bin][xp_bin],2)
                    cnt += 1
                if avg == 0: 
                    Avg[pT_bin][as_bin].append(-1)
                    AvgE[pT_bin][as_bin].append(0)
                    continue
                avg /= wei
                unc = 1./np.sqrt(wei)
                Avg[pT_bin][as_bin].append(avg)
                AvgE[pT_bin][as_bin].append(unc)
            AvgX[pT_bin][as_bin]=(X[0][pT_bin][as_bin])


plt.figure(1,figsize=(8,6))
plt.axvline(x=0.3,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.ylim([0.8,2.1])
plt.xlim([0.2,0.75])
plt.xticks([0.2,0.3,0.4,0.5,0.6,0.7],fontsize=16)
plt.yticks([1,1.5,2],fontsize=16)
plt.xlabel(r"$x'$",fontsize=22)
plt.ylabel('Data/Simulation',fontsize=22)




fig = plt.figure(figsize=(12, 9))
fig.subplots_adjust(hspace=0,wspace=0.07)

ax_loL= fig.add_subplot(3,2,5)  
ax_loR= fig.add_subplot(3,2,6,sharey=ax_loL)  
ax_miL= fig.add_subplot(3,2,3,sharex=ax_loL)  
ax_miR= fig.add_subplot(3,2,4,sharex=ax_loR,sharey=ax_miL)  
ax_upL= fig.add_subplot(3,2,1,sharex=ax_loL)  
ax_upR= fig.add_subplot(3,2,2,sharex=ax_loR,sharey=ax_upL)  

ax_loL.tick_params(direction='in')
ax_loR.tick_params(labelleft=False,direction='in')
ax_miL.tick_params(labelbottom=False,direction='in')
ax_miR.tick_params(labelbottom=False,labelleft=False,direction='in')
ax_upL.tick_params(labelbottom=False,direction='in')
ax_upR.tick_params(labelbottom=False,labelleft=False,direction='in')

fig.text(0.003, 0.5, 'Data/Simulation', va='center',rotation='vertical',fontsize=29)


plt.sca(ax_loL)
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.ylim([0.1,2.5])
plt.xlim([0.2,0.75])
plt.xticks([0.2,0.3,0.4,0.5,0.6,0.7],fontsize=16)
plt.yticks([0.5,1,1.5,2],fontsize=16)
plt.xlabel(r"$x'$",fontsize=22)

plt.sca(ax_miL)
plt.axvline(x=0.3,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.ylim([0.1,2.5])
plt.yticks([0.5,1,1.5,2],fontsize=16)

plt.sca(ax_upL)
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.axvline(x=0.3,linestyle='--',linewidth=2,color='black')
plt.ylim([0.1,2.5])
plt.yticks([0.5,1,1.5,2],fontsize=16)

plt.sca(ax_loR)
plt.xlim([0.2,0.75])
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.xticks([0.2,0.3,0.4,0.5,0.6,0.7],fontsize=16)
plt.xlabel(r"$x'$",fontsize=22)

plt.sca(ax_miR)
plt.axvline(x=0.3,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')

plt.sca(ax_upR)
plt.axvline(x=0.3,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')




files = ['ratio_10pt2.txt','ratio_10pt4.txt','ratio_10pt6.txt']
mark = ['o','s','*']
col = ['blue','red','green']

DataSimX    = [[] for j in range(len(files))] # array for each beam energy
DataSim     = [[] for j in range(len(files))]
DataSimErr  = [[] for j in range(len(files))]
Avg_DataSimX    = [[] for i in range(2) ]
Avg_DataSim     = [[] for i in range(2) ]
Avg_DataSimErr  = [[] for i in range(2) ]

for i in range(len(files)):
    DataSimX    [i] = [[] for j in range(2) ] # array for each pT setting
    DataSim     [i] = [[] for j in range(2) ]
    DataSimErr  [i] = [[] for j in range(2) ]
    for j in range(2):
        Avg_DataSimX    [j] = [ [] for k in range(3) ]
        Avg_DataSim     [j] = [ [] for k in range(3) ]
        Avg_DataSimErr  [j] = [ [] for k in range(3) ]
        DataSimX    [i][j] = [ [] for k in range(3) ] # array for each as setting
        DataSim     [i][j] = [ [] for k in range(3) ] # array for each as setting
        DataSimErr  [i][j] = [ [] for k in range(3) ] # array for each as setting

    with open(files[i],"r") as f:
        Q2_bin = -1
        pT_bin = -1
        as_bin = -1
        for line in f:
            if 'TCanvas' in line: continue
            if 'Files used' in line:
                parse = line.strip().split()
                check = float(parse[2].split("/")[4])
                against = 10+float(files[i].split("pt")[1].split(".")[0])/10
                if check != against: exit(-1)
                continue
            if 'Q2' in line:
                parse = line.strip().split(":")[1].split()
                Q2_bin = int(parse[0])
                pT_bin = int(parse[1])
                as_bin = int(parse[2])
                continue

            Read(line,DataSimX[i][pT_bin][as_bin], DataSim[i][pT_bin][as_bin], DataSimErr[i][pT_bin][as_bin] )

    
    for j in range(2):
        for k in range(3):
            
            #DataSimX[i][j][k] = np.asarray(DataSimX[i][j][k]) + (i-1)/200
            DataSimX[i][j][k] = np.asarray(DataSimX[i][j][k])
            '''
            if j == 0 and k == 0:

                plt.sca(ax_upL)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)
                plt.figure(1)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 0 and k == 1:
                plt.sca(ax_miL)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 0 and k == 2:
                plt.sca(ax_loL)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 0:
                plt.sca(ax_upR)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 1:
                plt.sca(ax_miR)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 2:
                plt.sca(ax_loR)
                plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)
            '''



WeightedAverage( Avg_DataSimX, Avg_DataSim, Avg_DataSimErr, DataSimX, DataSim, DataSimErr )


plt.sca(ax_upL)
plt.errorbar( Avg_DataSimX[0][0], Avg_DataSim[0][0] ,yerr=Avg_DataSimErr[0][0] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
plt.figure(1)
plt.errorbar( Avg_DataSimX[0][0], Avg_DataSim[0][0] ,yerr=Avg_DataSimErr[0][0] ,color='blue',marker='o',linestyle='none',markersize=12,linewidth=3)
plt.sca(ax_miL)
plt.errorbar( Avg_DataSimX[0][1], Avg_DataSim[0][1] ,yerr=Avg_DataSimErr[0][1] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
plt.sca(ax_loL)
plt.errorbar( Avg_DataSimX[0][2], Avg_DataSim[0][2] ,yerr=Avg_DataSimErr[0][2] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
plt.sca(ax_upR)
plt.errorbar( Avg_DataSimX[1][0], Avg_DataSim[1][0] ,yerr=Avg_DataSimErr[1][0] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
plt.sca(ax_miR)
plt.errorbar( Avg_DataSimX[1][1], Avg_DataSim[1][1] ,yerr=Avg_DataSimErr[1][1] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
plt.sca(ax_loR)
plt.errorbar( Avg_DataSimX[1][2], Avg_DataSim[1][2] ,yerr=Avg_DataSimErr[1][2] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)

plt.figure(1)
#l = plt.legend(numpoints=1,loc='best',ncol=1,framealpha=1,fontsize=15)
#l.set_zorder(99)



plt.figure(1)
ax=plt.gca()
textstr = '\n'.join((
    r'$p_e > %1.0f$ GeV/c' % (3, ),
    r'$Q^2 > %1.0f$ GeV$^2$/c$^4$' % (2, ),
    r'$W > %1.0f$ GeV/c$^2$' % (2, ),
    r'$E_\mathrm{dep} > %1.0f$ MeVee' % (10, ),
    r'$\cos\theta_{nq} < %1.1f$' % (-0.8, ),
    r'$p_n > %1.2f$ GeV/c' % (0.25, ), 
    r"$W' > %1.1f$ GeV/c$^2$" % (1.8, ),
    r'$p_T < %1.1f$ GeV/c' % (0.1, ),
    r'$ %1.1f < \alpha_S < %1.1f$' % (1.3,1.4, ),
    ))
props = dict(boxstyle='round', facecolor='white', alpha=1,zorder=99)
#plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=16,verticalalignment='top', bbox=props)
#plt.savefig('ratio_individual_firstbin.pdf',bbox_inches='tight')
plt.savefig('ratio_combo_firstbin.pdf',bbox_inches='tight')

#print("# x' \t As \t Pt \t Rat \t P2P ")
#for i in range(len(Avg_DataSimX[0][0])):
#    print(Avg_DataSimX[0][0][i],"\t",1.35,"\t",0.05,"\t",np.round(Avg_DataSim[0][0][i],3),"\t",np.round(Avg_DataSimErr[0][0][i],3))


#plt.savefig('ratio_individual.pdf',bbox_inches='tight')
#plt.savefig('ratio_combo_all.pdf',bbox_inches='tight')
plt.show()




