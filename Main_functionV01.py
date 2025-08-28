#Main function
from os import TMP_MAX
import numpy as np
import cmath
from matplotlib.patches import Ellipse

from generalCoefs import inputCoefs
from FractureEnergy import FractureEnergy
from ThermalExpansion import ThermalExpansion
from ETAMUClass import EtaMu
from VectorMatrixform import VectorMatrixForm
from abmnpq_Class import abmnpq
from TemperatureInclusions import TemperatureInclusions
from Matrix_Class import Matrix
from Plot_Class import plots


'initialize coeficients'

general = inputCoefs() 
fractureE = FractureEnergy() 
TE = ThermalExpansion(fractureE.TS)
EtMu = EtaMu(general)
VCM = VectorMatrixForm()
ab = abmnpq(general,EtMu)
Temp_In = TemperatureInclusions(TE , general , ab , VCM)
pt = plots(general)


#Test EtaMu:
"""print("EtaClosePQ = ",EtMu.EtaClosePQ(1,1,1,0))
print ("EtafarPQ = " , EtMu.EtafarPQ(1,1,1,0))
print("Etaint = "  , EtMu.EtaIntgPQ(1,1,1,0))
print("MuClosePQ = ",EtMu.MuClosePQ(1,1,1,0))
print ("MufarPQ = " , EtMu.MufarPQ(1,1,1,0))
print("Muint = "  , EtMu.MuIntgPQ(1,1,1,0))"""

#Test abmnpq
"""
i = 1
j =1
p=0
q=1
print("smnpq : " , ab.smnpq(i,j,p,q)) #OK
print("bmnpqa : " , ab.bmnpqa(i,j,p,q)) # OK
print("bmnpqb : " , ab.bmnpqb(i,j,p,q)) #OK
print("anpq = " , ab.anpq(i,p,q)) #OK if all amp equal to 1
print("bnpq = " , ab.bnpq(i,p,q)) # ok
print("bnpqminus = " , ab.bnpqminus(i,p,q)) # ok It deppends on a[i,j]
print ("dminus=" , ab.dminus(i,q)) #ok and it deppends on d[i,q]"""


#Temperature inclusions...


#print("M1Tem = " ,Temp_In.M1Tem) #Ok
#print("RealMTem array = " , Temp_In.RealMTem_array)
#print(Temp_In.realSTem)
#print("ImBCTem = " , Temp_In.ImBCTem)
#print("STem = " , Temp_In.STem)
#print("Snq = " , Temp_In.Snq[0,0])
#print("rmnpq = " , Temp_In.Rmnpq)  #ok
#print("f1 = " , Temp_In.f1)    #ok
#print("fnq1 = ",  Temp_In.fnq1[1][0])  #ok
#print("Fn = " , Temp_In.Fn(2,0))   #ok
#print("fnq = " , Temp_In.fnq(2,0)) #ok
#print("Gnq2 = " , Temp_In.Gnq2(1,0))

"""def Gnq (m,q):
    x = complex(np.real(Temp_In.Gnq4(m,q)) , np.imag(Temp_In.Gnq2(m,q)))
    return x"""

#print(Temp_In.f1)
#Gnq2array = VCM.VecCreate(Temp_In.Gnq2 , general.n , general.ntot)
#print("Gnq2array = " , Gnq2array)

#Gnq4array = VCM.VecCreate(Temp_In.Gnq4 , general.n , general.ntot)
#print("Gnq4array = " , Gnq4array)

#Gnqarray = VCM.VecCreate(Gnq , general.n , general.ntot)
#print("Gnq = " , Gnqarray)

#print("OmegaFar = " ,Temp_In.OmegaFar(1))
#print("OmegaSQ = " ,Temp_In.OmegaSQ(1))

#print("Xi = " , Temp_In.Xi)
#print("OmegaGQ = " ,Temp_In.OmegaGQ(1))
#print("Tfar = " , Temp_In.Tfar)
#print("Tmatrix = " , Temp_In.Tmatrix)
#print("TQ = " , Temp_In.TQ)
#print("Tfield total = " , Temp_In.Tfieldtotal)
#print("bcTem1nq = " , Temp_In.bcTem1nq(2,0))
print("bcTem2nq = " , Temp_In.bcTem2nq(2,0))