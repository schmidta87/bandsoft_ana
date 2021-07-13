import numpy as np
import matplotlib.pyplot as plt

barID = []
chi2 = []
colors=[]
with open("chi2_vals.txt","r") as f:
    for line in f:
        parse = line.strip().split()
        layer = int(parse[0])
        sector = int(parse[1])
        component = int(parse[2])
        barID.append(sector*100+layer*10+component)
        chi2.append( float(parse[3]) )
        if(float(parse[3]) > 6 or (float(parse[3]) != float(parse[3])) ):
            colors.append('red')
            print(layer,sector,component)
        else:
            colors.append('blue')

plt.scatter(barID,chi2,color=colors,marker='o')
plt.grid(True)
plt.xlim([100,570])
plt.ylim([0,20])
plt.savefig("chi_perbar.pdf")
plt.show()
