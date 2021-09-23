import numpy as np
import matplotlib.pyplot as plt



files = ["../bin/xdist_10pt2_newmixing.txt","../bin/xdist_10pt4_testnewbins.txt","../bin/xdist_10pt6_testnewbins.txt"]

chi2 = [ {} for i in files ]
colors = [ {} for i in files ]
marks = ['o','*','s']
maxes = []
plot_me = {}

for fi in range(len(files)):
    this_max = 0
    bad_bars = 0
    with open(files[fi],"r") as f:
        for line in f:
            if 'TCanvas' in line: continue
            if 'bar details' in line: continue
            parse = line.strip().split()

            if(len(parse)!=6): 
                print("error")
                print(line)
                exit(-1)
    
            layer = int(parse[2])
            sector = int(parse[3])
            comp = int(parse[4])
            ID = sector*100 + layer*10 + comp
            this_chi2 = float(parse[5])
            if this_chi2 > this_max: this_max = this_chi2 
            chi2[fi][ID] = this_chi2


            if this_chi2 > 2.0:
                print(fi,ID)
                bad_bars += 1
                plot_me[int(ID)] = True



           # if( this_chi2 < 2 ):
           #     plt.figure(1)
           #     plt.scatter( ID , this_chi2 , color = 'blue', marker = marks[fi] )
           # else:
           #     plt.figure(1)
           #     plt.scatter( ID , this_chi2 , color = 'red', marker = marks[fi] )
    print(bad_bars)
    maxes.append(this_max)

already_plotted = {}
for key in plot_me:
    already_plotted[key] = False
print(len(chi2[0]))
cnt = 0
for i in range(len(files)):
    for bar in chi2[i].keys():

        try:
            if plot_me[bar] == True and already_plotted[bar] == False:
                    cnt +=1
                    already_plotted[bar] = True
                    ymi = min([chi2[0][bar],chi2[1][bar],chi2[2][bar]])
                    yma = max([chi2[0][bar],chi2[1][bar],chi2[2][bar]])
                    print(bar,ymi,yma)

                    col = ''
                    plt.figure(1)
                    if chi2[0][bar] > 2.0: col = 'red'
                    else: col = 'blue'
                    plt.scatter( bar , chi2[0][bar] , color = col, marker = marks[0] )
                    if chi2[1][bar] > 2.0: col = 'red'
                    else: col = 'blue'
                    plt.scatter( bar , chi2[1][bar] , color = col, marker = marks[1] )
                    if chi2[2][bar] > 2.0: col = 'red'
                    else: col = 'blue'
                    plt.plot([bar,bar],[ymi,yma],color='black',linewidth=2,linestyle='--',zorder=-1)
                    plt.scatter( bar , chi2[2][bar] , color = col, marker = marks[2] ,linewidth=3)
        except KeyError:
            continue
print(cnt)
#print(cnt/len(files))

plt.figure(1)
plt.axhline(y=2.0,linestyle='--',color='black')
plt.legend(numpoints=1,loc='best')
plt.grid(True)
plt.xlim([100,570])
plt.ylim([0,4])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig("chi_perbar.pdf")
#
#
#
#
#
#
plt.show()
