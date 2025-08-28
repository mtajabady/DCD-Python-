from cmath import exp
import cmath
import numpy as np
from generalCoefs import inputCoefs
from ThermalExpansion import ThermalExpansion
from FractureEnergy import FractureEnergy
from TemperatureInclusions import TemperatureInclusions 






class farfield:
    def __init__(self , TE : ThermalExpansion ,  general : inputCoefs , fractureE : FractureEnergy , Temp_In : TemperatureInclusions):
        self.general = general
        self.TE = TE
        self.aa = []
        self.fractureE = fractureE
        self.Temp_In = Temp_In
        for q in range (self.general.ntot):
            self.aa.append(cmath.exp(complex(0,-self.general.Theta[q]))*(self.general.dd(q)/(16*self.general.Mu0))*(self.general.s11+self.general.s22))
        self.bbm = []
        for q in range(self.general.ntot):
            self.bbm.append(self.aa[q]*self.general.v0(0)**(-2) + cmath.exp(complex(0,self.general.Theta[q]))*(self.general.dd(q)/(8*self.general.Mu0))*(self.general.s22-self.general.s12 + 2*complex(0,self.general.s12)))
            "faraa[q_]  := Module[{f}, f= Table[If[i ==1&&q<= ntot, Exp[-I Theta[q]](dd[q]/(16  Mu0))(s11+s22),0],{i,1,n}]]"            
        self.r5 = []
        for q in range(self.general.ntot):
            self.r5.append(self.faraa(q))
        self.faraaminus = []
        for q in range (self.general.ntot):
            self.faraaminus.append(self.r5[q])
        self.listaaminus =[]
        for q in range (self.general.ntot):
            self.listaaminus.append(self.listaa(q))
        self.bcr = np.real(self.bcc())
        self.bci = np.imag(self.bcc())
        self.bcbar = [self.bcr , self.bci]

        
        
        pass
    def faraa(self,q):
        f = []
        for i in range (1,self.general.n+1):
            if i==1 and q<= self.general.ntot :
                f.append((cmath.exp(complex(0,-self.general.Theta[q]))*(self.general.dd(q)/(16*self.general.Mu0))*(self.general.s11+self.general.s22)))
            else:
                f.append(0)
        return f
    def farbbminus(self,q):
        f = []
        f1 = (cmath.exp(complex(0,-self.general.Theta[q]))*(self.general.dd(q)/(16*self.general.Mu0))*(self.general.s11+self.general.s22))
        f2 = f1*self.general.v0(q)**(-2) + (cmath.exp(complex(0,self.general.Theta[q]))*(self.general.dd(q)/(8*self.general.Mu0))*(self.general.s22-self.general.s11+2*complex(0,self.general.s12)))
        for i in range (1,self.general.n+1):
            if i==1 and q<=self.general.ntot:
                f.append(f2)
            else:
                f.append(0)
        return f
    def farbb(self,q):
        f=[]
        for i in range (1,self.general.n+1):
            f.append(self.farbbminus(q)[i-1] + 2*i*cmath.sinh(2*self.general.Zeta0(q))*self.faraa(q)[i-1])
        return f
    def listaa(self,q):
        l=[]
        for i in range(1,self.general.n+1):
            l.append(self.faraa(q)[i-1])
        return l
    def listbbminus(self,q):
        b=[]
        for i in range (1, self.general.n + 1):
            b.append(self.farbbminus(q)[i-1])
        return b
    def listbb(self , q):
        b=[]
        for i in range(1,self.general.n+1):
            b.append(self.farbb(q)[i-1])
        return b
    def bcf1(self,i,j,q):
        f1 = -self.general.Chi0*self.listaa(q)[j-1] + (self.general.v0(q)**(2*j))*np.conj(self.listbbminus(q)[j-1])
        if self.fractureE.TS == 0 :
            f = f1
        else:
            f = f1 + self.Temp_In.bcTem1nq(j,q)
        if i==1:
            x=f
        else:
            x=0
        return x

    def bcf2(self,i,j,q):
        f1 = -self.general.Chi0*self.general.v0(q)**(2*j)*self.listaaminus[q][j-1] - np.conj(self.listbb(q)[j-1])
        if self.fractureE.TS == 0 :
            f = f1
        else:
            f = f1 + self.Temp_In.bcTem2nq(j,q)
        if i==2:
            x=f
        else:
            x=0
        return x
    def bcf3(self,i,j,q):
        if i==3:
            f = -self.listaa(q)[j-1] - self.general.v0(q)**(2*j)*np.conj(self.listbbminus(q)[j-1])
        else:
            f = 0
        return f
    def bcf4(self,i,j,q):
        if i==4:
            f = -self.general.v0(q)**(2*j)*self.listaaminus[q][j-1] - np.conj(self.listbb(q)[j-1])
        else:
            f = 0
        return f
    def bcf(self,i,j,q):
        return self.bcf1(i,j,q) + self.bcf2(i,j,q) + self.bcf3(i,j,q) + self.bcf4(i,j,q)
    def bbc(self,j,q):
        b = []
        for i in range(1,self.general.n+1):
            b.append(self.bcf(j,i,q))
        return b
    def bc(self,q):
        bc=[]
        for i in range(1,4):
            bc.append(self.bbc(i,q))
        return bc
    def bbctest(self,j,q):
        b = []
        for i in range(1,5):
            b.append(self.bcf(j,i,q))
        return b
    def bcc(self):
        for q in range (self.general.ntot):
            for i in range(1,self.general.n + 1):
                f = self.bbctest(i,q)
        return f
    




        
