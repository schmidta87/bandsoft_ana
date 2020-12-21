import numpy as np
import matplotlib.pyplot as plt

files = ["tdc_shifts.txt","fadc_shifts.txt"]

i = 0
for fi in files:
    IDs = []
    Shifts = []
    Sigmas = []
    with open("../bin/"+fi,"r") as f:
        for line in f:
            parse = line.strip().split()
            sector = int(parse[0])
            layer = int(parse[1])
            component = int(parse[2])

            shift = float(parse[5])
            sigma = float(parse[6])

            IDs.append( layer*10 + sector*100 + component )
            Shifts.append( shift )
            Sigmas.append( sigma )

    plt.figure(i)
    plt.errorbar( IDs, Shifts , yerr = Sigmas , linestyle='none' , marker='o' )
    plt.savefig('shifts_'+str(i)+'.pdf',bbox_inches='tight')
    i+=1

plt.show()
