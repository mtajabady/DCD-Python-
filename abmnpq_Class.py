import cmath
from generalCoefs import inputCoefs
from ETAMUClass import EtaMu 
import numpy as np


class abmnpq:
    def __init__(self, general : inputCoefs, EM : EtaMu) :
        self.general = general
        self.EtaMu = EM
        self.amp = np.ones(general.n*general.ntot).reshape(general.n,general.ntot)
        self.bmp = np.ones(general.n*general.ntot).reshape(general.n,general.ntot)
        self.cmp = np.ones(general.n*general.ntot).reshape(general.n,general.ntot)
        self.dmp = np.ones(general.n*general.ntot).reshape(general.n,general.ntot)

        pass       
    def smnpq(self,m,i,p,q):
        x = self.EtaMu.EtaIntgPQ(m , abs(i) , p ,q)*cmath.exp(complex(0,-self.general.fTheta(p,q)))
        return x
    def sum_c(self,l,i,m,p,q):                    #   calculation of M_nml for k = 0 to k = l
        sum01 = 0
        for k in range(l):
            x= (2*k + abs(i))*self.EtaMu.EtaIntgPQ(m,2*k+abs(i) , p,q)
            sum01 = sum01 + x
        return (sum01)

    def bmnpqa (self,m,i,p,q):   
        if m ==1:
            e_if = self.EtaMu.MuClosePQ(m-1,abs(i),p,q)
        else:
            e_if = self.EtaMu.MuIntgPQ(m-1,abs(i),p,q)
        if p==q:
            final=0
        else:
            a = (-self.general.dd(p)/2)* (m/(m+1)) * ((self.general.v0(p) - 1/self.general.v0(p))**2) * cmath.exp(complex(0,self.general.fTheta(p,q)))*self.EtaMu.MuIntgPQ(m+1,abs(i),p,q)
            b = (abs(i)*cmath.exp(complex(0,-self.general.fTheta(p,q)))*(-self.general.v0(q)**(-2)+1)-m*cmath.exp(complex(0,self.general.fTheta(p,q)))*(-(self.general.v0(p)**(-2))+1))*self.EtaMu.EtaIntgPQ(m,abs(i),p,q) #OK
            c = ((self.general.v0(q)-1/self.general.v0(q))**2)*cmath.exp(complex(0,-self.general.fTheta(p,q)))*self.sum_c(self.general.nmax,i,m,p,q)
            d = (cmath.exp(complex(0,self.general.Theta[q]))*np.conj(self.EtaMu.zpq(p,q)) - self.EtaMu.zpq(p,q)*cmath.exp(complex(0,-self.general.Theta[q])))*cmath.exp(complex(0,self.general.Theta[p]))*self.EtaMu.MuIntgPQ(m,abs(i),p,q)
            e = ((self.EtaMu.MuIntgPQ(m+1 , abs(i) , p,q) + e_if)*self.general.dd(p)/2 -self.EtaMu.EtaIntgPQ(m,abs(i),p,q))*(cmath.exp(complex(0 , -self.general.fTheta(p,q)))-cmath.exp(complex(0 , self.general.fTheta(p,q))))
            final = a+b+c+d+e
        return final
    def bmnpqb (self,m,i,p,q):   
        f = self.EtaMu.EtaIntgPQ(m,abs(i),p,q)*cmath.exp(complex(0,self.general.fTheta(p,q)))
        return f
    def anpq(self,i,p,q):
        sum01 = 0
        for m in range(self.general.n):
            x= self.amp[m,p]*self.smnpq(m,i,p,q)
            sum01 = sum01 + x
        return (sum01)
    def sumbnpq(self,i,p,q):
        sum01 = 0
        for m in range(self.general.n):
            x= self.bmp[m,p]*self.bmnpqb(m,i,p,q)+self.amp[m,p]*self.bmnpqa(m,i,p,q)
            sum01 = sum01 + x
        return (sum01)
    def bnpq(self,i,p,q):
        if p==q:
            f=0
        else:
            f = self.sumbnpq(i,p,q)
        return f
    def bnpqminus(self,i,p,q):
        f = self.bnpq(i,p,q)-2*i*cmath.sinh(2*self.general.Zeta0(q))*self.anpq(i,p,q)
        return f
    def dminus(self,i,q):
        f = self.dmp[i,q]-2*i*cmath.sinh(2*self.general.Zeta0(q))*self.cmp[i,q]
        return f
             
