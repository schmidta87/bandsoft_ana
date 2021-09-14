import numpy as np
import matplotlib.pyplot as plt

def Read(line,X,Y,E):
    parse = line.strip().split()
    x = float(parse[1])
    y = float(parse[2])
    e = float(parse[3])
    X.append(x)
    Y.append(y)
    E.append(e)

fig, (aU, aL) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(10,8))
fig.subplots_adjust(hspace=0,wspace=0)

aU.tick_params(labelbottom=False,direction='in')
aL.tick_params(direction='in')


#plt.figure(1,figsize=(9,6))
plt.sca(aL)
plt.axhline(y=1,linestyle='--',linewidth=3,color='black')
plt.ylim([0.5,1.5])
plt.xlim([1.30,1.6])
plt.xticks([1.3,1.4,1.5,1.6],fontsize=14)
plt.yticks([0.7,1,1.3],fontsize=14)
plt.ylabel('Data/Simulation',fontsize=18,labelpad=12)
plt.xlabel(r'$\alpha_S$',fontsize=18)

#plt.figure(2,figsize=(9,6))
plt.sca(aU)
plt.xlim([1.30,1.6])
plt.ylim([0,1800])
plt.xticks([1.3,1.4,1.5,1.6],fontsize=14)
plt.yticks([0,500,1000,1500],fontsize=14)
plt.ylabel('Counts [a.u.]',fontsize=18,labelpad=12)
plt.xlabel(r'$\alpha_S$',fontsize=18)

files = ['as_10pt2.txt','as_10pt4.txt','as_10pt6.txt']
mark = ['o','s','*']
col = ['blue','red','green']

DataSimX    = [[] for j in range(len(files))]
DataSim     = [[] for j in range(len(files))]
DataSimErr  = [[] for j in range(len(files))]
DataX       = [[] for j in range(len(files))]
Data        = [[] for j in range(len(files))]
DataErr     = [[] for j in range(len(files))]
SimX        = [[] for j in range(len(files))]
Sim         = [[] for j in range(len(files))]
SimErr      = [[] for j in range(len(files))]

for i in range(len(files)):
    DataSimX    [i] = [[] for j in range(3) ]
    DataSim     [i] = [[] for j in range(3) ]
    DataSimErr  [i] = [[] for j in range(3) ]
    DataX       [i] = [[] for j in range(3) ]
    Data        [i] = [[] for j in range(3) ]
    DataErr     [i] = [[] for j in range(3) ]
    SimX        [i] = [[] for j in range(3) ]
    Sim         [i] = [[] for j in range(3) ]
    SimErr      [i] = [[] for j in range(3) ]

    
    j = -1
    with open(files[i],"r") as f:
        for line in f:
            if 'Files used' in line:
                parse = line.strip().split()
                check = float(parse[2].split("/")[4])
                against = 10+float(files[i].split("pt")[1].split(".")[0])/10
                if check != against: exit(-1)
                continue
            if 'pT' in line:
                j += 1

            if 'Data/Sim:' in line:
                Read(line,DataSimX[i][j],DataSim[i][j],DataSimErr[i][j])
            elif 'Data:' in line:
                Read(line,DataX[i][j],Data[i][j],DataErr[i][j])
            elif 'Sim:' in line:
                Read(line,SimX[i][j],Sim[i][j],SimErr[i][j])
            else: continue

        plt.sca(aL)
        DataSimX[i][1] = np.asarray(DataSimX[i][1]) + (i-1)/300
        plt.errorbar(DataSimX[i][1],DataSim[i][1],yerr=DataSimErr[i][1],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV')

        plt.sca(aU)
        DataX[i][1] = np.asarray(DataX[i][1]) + (i-1)/300
        plt.errorbar(DataX[i][1],Data[i][1],yerr=DataErr[i][1],linestyle='none',marker=mark[i],color=col[i],label=str(check)+' GeV')

        SimX[i][1] = np.asarray(SimX[i][1])+0.0125
        Sim[i][1] = np.asarray(Sim[i][1])
        SimErr[i][1] = np.asarray(SimErr[i][1])
        plt.step(SimX[i][1],Sim[i][1],color=col[i])
        #plt.step(SimX[1],Sim[1]+SimErr[1],color=col[i],alpha=0.5)
        #plt.step(SimX[1],Sim[1]-SimErr[1],color=col[i],alpha=0.5)


Avg_DataSimX    = []
Avg_DataSim     = []
Avg_DataSimErr  = []
Avg_DataX       = []
Avg_Data        = []
Avg_DataErr     = []
Avg_SimX        = []
Avg_Sim         = []
Avg_SimErr      = []
def WeightedAverage( AvgX, Avg, AvgE, X, Y, E ):
    for i in range(len(X[0][1])): # for each bin in X:

        avg = 0
        wei = 0
        cnt = 0
        for j in range(3):
            if E[j][1][i] == 0: continue
            avg += Y[j][1][i]/pow(E[j][1][i],2)
            wei += 1./pow(E[j][1][i],2)
            cnt += 1
        if avg == 0: 
            Avg.append(-1)
            AvgE.append(0)
            AvgX.append(X[j][1][i])
            continue
        avg /= wei
        unc = 1./np.sqrt(wei)
        Avg.append(avg)
        AvgE.append(unc)
        AvgX.append(X[j][1][i])


WeightedAverage( Avg_DataSimX, Avg_DataSim, Avg_DataSimErr, DataSimX, DataSim, DataSimErr )

WeightedAverage( Avg_DataX, Avg_Data, Avg_DataErr, DataX, Data,       DataErr )

WeightedAverage( Avg_SimX, Avg_Sim,  Avg_SimErr, SimX, Sim,       SimErr )




plt.sca(aL)
#plt.scatter( Avg_DataSimX, Avg_DataSim ,color='blue',marker='x',linewidth=4)
#plt.errorbar( Avg_DataSimX, Avg_DataSim ,yerr=Avg_DataSimErr ,color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)

plt.sca(aU)
#plt.scatter( Avg_DataX, Avg_Data ,color='blue',marker='x',linewidth=4)
#plt.errorbar( Avg_DataX, Avg_Data ,yerr=Avg_DataErr, color='blue',marker='x',linestyle='none',markersize=12,linewidth=3)
#plt.step( Avg_SimX , Avg_Sim , color='black')

plt.legend(numpoints=1,loc=1,fontsize=12)


plt.savefig('data_as_individual.pdf',bbox_inches='tight')
#plt.savefig('data_as_combo.pdf',bbox_inches='tight')
plt.show()
