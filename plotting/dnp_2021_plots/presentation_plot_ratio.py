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
        print(pT_bin,len(X[0]))
        for as_bin in range(len(X[0][pT_bin])):
            print(as_bin,len(X[0][pT_bin]))
            for xp_bin in range(len(X[0][pT_bin][as_bin])):
                print(xp_bin,len(X[0][pT_bin][as_bin]))
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
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.ylim([0.1,2.5])
plt.yticks([0.5,1,1.5,2],fontsize=16)

plt.sca(ax_upL)
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.ylim([0.1,2.5])
plt.yticks([0.5,1,1.5,2],fontsize=16)

plt.sca(ax_loR)
plt.xlim([0.2,0.75])
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.xticks([0.2,0.3,0.4,0.5,0.6,0.7],fontsize=16)
plt.xlabel(r"$x'$",fontsize=22)

plt.sca(ax_miR)
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')

plt.sca(ax_upR)
plt.axvline(x=0.4,linestyle='--',linewidth=2,color='black')
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
            
            DataSimX[i][j][k] = np.asarray(DataSimX[i][j][k]) + (i-1)/200
            if j == 0 and k == 0:
                plt.sca(ax_upL)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 0 and k == 1:
                plt.sca(ax_miL)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 0 and k == 2:
                plt.sca(ax_loL)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 0:
                plt.sca(ax_upR)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 1:
                plt.sca(ax_miR)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

            if j == 1 and k == 2:
                plt.sca(ax_loR)
                #plt.errorbar(DataSimX[i][j][k],DataSim[i][j][k],yerr=DataSimErr[i][j][k],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)


            #if 'Data/Sim:' in line:
                #Read(line,DataSimX[i][j],DataSim[i][j],DataSimErr[i][j])
            #elif 'Data:' in line:
                #Read(line,DataX[i][j],Data[i][j],DataErr[i][j])
            #elif 'Sim:' in line:
                #Read(line,SimX[i][j],Sim[i][j],SimErr[i][j])
            #else: continue
        
        #plt.sca(aL)
        #DataSimX[i][0] = np.asarray(DataSimX[i][0]) + (i-1)/1000
        #plt.errorbar(DataSimX[i][0],DataSim[i][0],yerr=DataSimErr[i][0],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=7)

        #plt.sca(aU)
        #DataX[i][0] = np.asarray(DataX[i][0]) + (i-1)/1000
        #plt.errorbar(DataX[i][0],Data[i][0],yerr=DataErr[i][0],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV',markersize=9)

        #SimX[i][0] = np.asarray(SimX[i][0])+0.01/2
        #Sim[i][0] = np.asarray(Sim[i][0])
        #SimErr[i][0] = np.asarray(SimErr[i][0])
        #plt.step(SimX[i][0],Sim[i][0],color=col[i],linewidth=2,zorder=-1)
        #plt.step(SimX[1],Sim[1]+SimErr[1],color=col[i],alpha=0.5)
        #plt.step(SimX[1],Sim[1]-SimErr[1],color=col[i],alpha=0.5)


WeightedAverage( Avg_DataSimX, Avg_DataSim, Avg_DataSimErr, DataSimX, DataSim, DataSimErr )


plt.sca(ax_upL)
plt.errorbar( Avg_DataSimX[0][0], Avg_DataSim[0][0] ,yerr=Avg_DataSimErr[0][0] ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
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


l = plt.legend(numpoints=1,loc=2,ncol=3,framealpha=1)
l.set_zorder(99)




#plt.savefig('ratio_individual.pdf',bbox_inches='tight')
plt.savefig('ratio_combo.pdf',bbox_inches='tight')
plt.show()




exit(-1)
'''


#WeightedAverage( Avg_DataSimX, Avg_DataSim, Avg_DataSimErr, DataSimX, DataSim, DataSimErr )

#WeightedAverage( Avg_DataX, Avg_Data, Avg_DataErr, DataX, Data,       DataErr )

#WeightedAverage( Avg_SimX, Avg_Sim,  Avg_SimErr, SimX, Sim,       SimErr )


cuts_str = '\n'.join((
    r'$\mu=%.2f$' % (2, ),
    r'$\mathrm{median}=%.2f$' % (2, ),
    r'$\sigma=%.2f$' % (2, )))

plt.sca(aL)
#plt.scatter( Avg_DataSimX, Avg_DataSim ,color='blue',marker='x',linewidth=4)
#plt.errorbar( Avg_DataSimX, Avg_DataSim ,yerr=Avg_DataSimErr ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)

plt.sca(aU)
#plt.text(0.05, 0.95, cuts_str, transform=ax.transAxes, fontsize=14,
#        verticalalignment='top', bbox=props)
#plt.scatter( Avg_DataX, Avg_Data ,color='blue',marker='x',linewidth=4)
#plt.errorbar( Avg_DataX, Avg_Data ,yerr=Avg_DataErr, color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
#plt.step( Avg_SimX , Avg_Sim , color='black')

plt.legend(numpoints=1,loc=2,fontsize=16)


plt.savefig('data_pT_individual.pdf',bbox_inches='tight')
#plt.savefig('data_pT_combo.pdf',bbox_inches='tight')
plt.show()
'''
