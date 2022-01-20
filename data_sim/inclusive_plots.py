import numpy as np
import matplotlib.pyplot as plt
import os

directory = '/Volumes/BAND-Data/v3.1/10.2/final/inclusive/fiducial_comparisons'

theta_bins = 7
min_theta = 8.5
step_theta = 3
mom_bins = 9
min_mom = 3
step_mom = 0.5

Phi = [ [ {} for j in range(mom_bins)] for i in range(theta_bins) ]

for fname in os.listdir(directory):
    if not '.txt' in fname: continue
    
    filename = directory+"/"+fname
    print('working on: ',filename)
    run = int(filename.split("/")[-1].split("_")[0])
    with open( filename, "r" ) as f:
        for line in f:
            parse = line.strip().split()

            if parse[0] == 'Phi':
                this_theta      = float(parse[2])
                this_theta_bin  = int((this_theta - min_theta)/step_theta)

                this_mom        = float(parse[8])
                this_mom_bin    = int((this_mom - min_mom)/step_mom)
                
                test_stat       = float(parse[15])
                dof             = int(parse[18])

                Phi[this_theta_bin][this_mom_bin][run] = [test_stat,dof,line.strip()] 
            else:
                continue

for i in range(theta_bins):
    for j in range(mom_bins):
        X = []
        Y1 = []
        Y2 = []
        info = ""
        for run in list(Phi[i][j].keys()):
            X.append(run)
            Y1.append(Phi[i][j][run][0])
            Y2.append(Phi[i][j][run][1])
            info = Phi[i][j][run][2]
        if( X or Y1 or Y2 ):
            plt.figure(i*100+j,figsize=(10,8))
            plt.scatter(X,Y1)
            
            parse = info.split()
            plt.xlabel('Run Number',fontsize=16)
            plt.ylabel('Test statistic',fontsize=16)
            title = ""
            save = ""
            for k in range(0,13):
                title += (parse[k] + " ")
            plt.title(title)
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)
            plt.axvline(x=6420,color='red',linestyle='--')

            plt.savefig(parse[0]+"_"+parse[2]+"_"+parse[6]+"_"+parse[8]+"_"+parse[12]+".pdf")

            plt.close(i*100+j)

