import sys
from F1F2 import f1, f2
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d
#import pystyle

me = 0.000511
mU = 0.9314941024
mN = 0.93892
mD = 2.01410178 * mU - me;
mbar = mD/2.

GeVfm = 0.1973;
alphaEM = 0.0072973525664;
cmSqGeVSq = GeVfm*GeVfm*1.E-26
nbGeVSq = cmSqGeVSq*1.E33

k = []
phiSq = []

with open("deuteron.csv") as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        split = row[0].split()
        k.append(float(split[0]))
        phiSq.append(float(split[3]))

k = np.array(k)*GeVfm
phiSq = np.array(phiSq)/(2*np.pi*GeVfm)**3

n = interp1d(k,phiSq,bounds_error = False,fill_value = 0)

def rho(alphaS,kt):
    kSq = (kt**2 + mN**2)*(alphaS*(2.-alphaS))**-1 - mN**2
    return np.sqrt(kSq+mN**2)*(2-alphaS)**-1*n(np.sqrt(kSq))
    
def F1p(x,qsq):
    wsq = qsq*(1./x-1.)+mN*mN
    return f1(1,1,qsq,wsq)

def F2p(x,qsq):
    wsq = qsq*(1./x-1.)+mN*mN
    return f2(1,1,qsq,wsq)

def x_QSq_tilde(xB,QSq,alphaS,kt):
    omega = QSq*(2.*mN*xB)**(-1)
    q = np.sqrt(QSq+omega**2)
    qplus = omega - q
    qminus = omega + q
    pplusS = alphaS*mbar
    pplusI = (2.-alphaS)*mbar
    pminusS = (mN**2+kt**2)*pplusS**(-1)
    pminusI = (mN**2+kt**2)*pplusI**(-1) # On the mass shell

    qminus_tilde = qminus + mD - pminusI - pminusS # Violation of minus-momentum

    QSq_tilde = -qplus*qminus_tilde
    x_tilde = QSq_tilde/(qplus*pminusI + qminus_tilde*pplusI)

    return x_tilde, QSq_tilde

def asym(xB,QSq,alphaS,kt,E=10.2):
    omega = QSq/(2.*mN*xB)
    y = omega/E
    if (y < 0) or (y > 1):
        return 0
    x_tilde, QSq_tilde = x_QSq_tilde(xB,QSq,alphaS,kt)
    F1 = F1p(x_tilde,QSq_tilde)
    F2 = F2p(x_tilde,QSq_tilde)
    if (F1 < 0) or (F2 < 0):
        return 0
    
    return (nbGeVSq*2.*alphaEM**2)*(xB*QSq**2)**(-1)*(x_tilde*y**2*F1 + (1.-y-(xB*y*mN)**2*QSq**(-1))*F2)*rho(alphaS,kt)

def full(xB,QSq,alphaS,kt,E=10.2,phi = np.pi):
    # phi is the initial proton and the initial electron, around the momentum transfer
    omega = QSq/(2.*mN*xB)
    q = np.sqrt(QSq+omega**2)
    qplus = omega - q
    qminus = omega + q
    pplusS = alphaS*mbar
    pplusI = (2.-alphaS)*mbar
    pminusS = (mN**2+kt**2)*pplusS**(-1)
    pminusI = (mN**2+kt**2)*pplusI**(-1) # On the mass shell

    qminus_tilde = qminus + mD - pminusI - pminusS # Violation of minus-momentum
    y = omega/E
    if (y < 0) or (y > 1):
        return 0
    
    Ek = E - omega

    pe_z = -(omega*E + QSq/2.)/q
    pk_z = -(omega*Ek - QSq/2.)/q

    pe_plus = E + pe_z
    pk_plus = Ek + pk_z
    pe_minus = E - pe_z
    pk_minus = Ek - pk_z

    pe_perp = np.sqrt(pe_plus*pe_minus)
    
    x_tilde, QSq_tilde = x_QSq_tilde(xB,QSq,alphaS,kt)
    F1 = F1p(x_tilde,QSq_tilde)
    F2 = F2p(x_tilde,QSq_tilde)
    if (F1 < 0) or (F2 < 0):
        return 0

    Lplus = pplusI + qplus/(2.*x_tilde)

    wWA = 4.*pe_minus*pk_minus*(-qplus**2/QSq_tilde*F1+Lplus**2*2*x_tilde/QSq_tilde*F2)

    wWB = 4.*(pe_perp**2*F1+(pe_perp*kt*np.cos(phi))**2*2*x_tilde*QSq_tilde*F2)

    wWC = QSq*(2.*F1 + kt**2*2*x_tilde/QSq_tilde*F2)

    wWD = 4.*(pe_minus + pk_minus)*Lplus*pe_perp*kt*np.cos(phi)*2*x_tilde/QSq_tilde*F2
    
    wW = wWA + wWB + wWC + wWD
    

    xDxB_jacobian = mN/mD # additional factor to translate differential cross section from dx_d -> dx_B
    alphaS_jacobian = alphaS # multiply by alphaS to translate differential cross section from dAlphaS/alphaS -> dAlphaS
    return nbGeVSq * (2.*mN)/mD * (alphaEM**2*y**2)/((2.-alphaS)*QSq**3)*wW*rho(alphaS,kt) \
            * xDxB_jacobian * alphaS_jacobian #units of nb*GeV^-4



# Open a list of kinematic points and give back the differential cross section for those points!
ofile_dat = open("data_bin_avgpoints_pointcs.txt", "w")
ofile_sim = open("sim_bin_avgpoints_pointcs.txt", "w")


def writeoutput(infile,outfile):
    outfile.write("# [binQ2] [binPt] [binXb] [binAs] [Q2] [Pt] [Xb] [As] [DifferntialCS (nb*GeV^-4)]\n")
    with open(infile,"r") as f:
        for line in f:
            if '#' in line: continue
            parse = line.strip().split()
            binQ2 = int(parse[0])
            binPt = int(parse[1])
            binXb = int(parse[2])
            binAs = int(parse[3])

            Q2 = float(parse[4])
            Pt = float(parse[5])
            Xb = float(parse[6])
            As = float(parse[7])
            
            if( Q2 == 0 or Pt == 0 or Xb == 0 or As == 0 ): continue

            print(binQ2,binPt,binXb,binAs,Q2,Pt,Xb,As,asym(Xb,Q2,As,Pt,10.2) )
            outfile.write( str(binQ2) + " " + str(binPt) + " " + str(binXb) + " " + str(binAs) + " " +\
                        str(Q2) + " " + str(Pt) + " " + str(Xb) + " " + str(As) + " " + str(asym(Xb,Q2,As,Pt,10.2)) + "\n")
        
writeoutput("data_bin_avgpoints.txt",ofile_dat)
writeoutput("sim_bin_avgpoints.txt",ofile_sim)

