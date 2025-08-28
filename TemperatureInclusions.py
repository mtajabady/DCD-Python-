
from cmath import exp
import cmath
import numpy as np
import sympy as sym
from generalCoefs import inputCoefs
from ThermalExpansion import ThermalExpansion
from abmnpq_Class import abmnpq
from VectorMatrixform import VectorMatrixForm
from Matrix_Class import Matrix


class TemperatureInclusions:
    def __init__(self , TE : ThermalExpansion ,  general : inputCoefs , ab: abmnpq , VCM: VectorMatrixForm ):
        self.general = general
        self.TE = TE
        self.ab = ab
        self.VCM = VCM
        self.t0 = 20
        self.gamma = 21.213203435596427-21.213203435596427j
        self.x = 1
        self.y = 1
        self.z = complex(self.x , self.y)
        self.Xi =[]
        for q in range(self.general.ntot):
            self.Xi.append(cmath.acosh((self.z - self.general.zcenter[q])/self.general.dd(q)))
        self.Omegafar0 = []
        for i in range(self.general.ntot):
            self.Omegafar0.append(self.t0 + 0*np.conj(self.gamma)*self.general.zcenter[i])
        self.maxn = self.general.n
        self.Ffar = VCM.VecCreate(self.Fn , general.n , general.ntot)
        self.BCTem = VCM.VecCreate(self.BCTem1 , general.n , general.ntot)
        self.realBCTem =  np.real(self.BCTem)
        self.ImBCTem = np.imag(self.BCTem)
        self.m01=[]
        self.m1 = VCM.MatrixCreate(self.m01,general.n,general.ntot)
        self.M1Tem = Matrix(self.general.n*self.general.ntot , self.general.n*self.general.ntot)
        for i in range (1,self.general.n*self.general.ntot+1 ):
            for j in range (1,self.general.n*self.general.ntot+1 ):
                m = self.m1[i-1,j-1,0]
                n = self.m1[i-1,j-1,1]
                p = self.m1[i-1,j-1,2]-1
                q = self.m1[i-1,j-1,3]-1
                Mr = np.real(self.M1Tem1(m,n,p,q))
                Mi = np.imag(self.M1Tem1(m,n,p,q))
                Mc = complex(Mr,Mi)
                Matrix.setElement(self.M1Tem,i,j,Mc)
        self.m02=[]
        self.m2 = VCM.MatrixCreate(self.m02 , general.n , general.ntot)
        self.M2Tem = Matrix(self.general.n*self.general.ntot , self.general.n*self.general.ntot)
        " this for loop changes M2Tem indices with function M1Tem by same order of (m,n,p,q) in m02"
        "Attention: m , n start from 1 but p,q start from 0"
        for i in range (1,self.general.n*self.general.ntot+1 ):
            for j in range (1,general.n*general.ntot+1 ):
                m = self.m2[i-1,j-1,0]
                n = self.m2[i-1,j-1,1]
                p = self.m2[i-1,j-1,2]-1
                q = self.m2[i-1,j-1,3]-1
                Mr = np.real(self.M2Tem1(m,n,p,q))
                Mi = np.imag(self.M2Tem1(m,n,p,q))
                Mc = complex(Mr,Mi)
                Matrix.setElement(self.M2Tem,i,j,Mc)
        
        self.RealMTem = Matrix(self.general.n*self.general.ntot , self.general.n*self.general.ntot)
        for i in range(general.n*general.ntot):
            for j in range(general.n*general.ntot):
                Matrix.setElement(self.RealMTem , i , j , (self.M1Tem[i-1][j-1]+self.M2Tem[i-1][j-1])) #indices in matrix creation start from 1 but in reading start from 0

        self.ImMTem = Matrix(self.general.n*self.general.ntot , self.general.n*self.general.ntot)
        for i in range(1,self.general.n*self.general.ntot+1):
            for j in range(1,self.general.n*self.general.ntot+1):
                Matrix.setElement(self.ImMTem , i , j , self.M1Tem[i-1][j-1] - self.M2Tem[i-1][j-1])   #indices in matrix creation start from 1 but in reading start from 0

        " change matrix[i][j] to array[i,j] to use linear solve"
        self.RealMT_array0 = []
        for i in range(self.general.n*self.general.ntot):
            for j in range(self.general.n*self.general.ntot):
                self.RealMT_array0.append(self.RealMTem[i][j])
        self.RealMTem_array = np.array(self.RealMT_array0).reshape(self.general.n*self.general.ntot ,self.general.n*self.general.ntot )
        self.realSTem = np.linalg.solve(self.RealMTem_array, self.realBCTem)

        self.ImMT_array0 = []
        for i in range(self.general.n*self.general.ntot):
            for j in range(self.general.n*self.general.ntot):
                self.ImMT_array0.append(self.ImMTem[i][j])
        self.ImMTem_array = np.array(self.ImMT_array0).reshape(self.general.n*self.general.ntot ,self.general.n*self.general.ntot )
        self.ImSTem = np.linalg.solve(self.ImMTem_array, self.realBCTem)

        self.STem = []
        for i in range(len(self.realSTem)):
            self.STem.append(complex(self.realSTem[i] , -self.ImSTem[i]))

        self.Snq1 = []
        self.nmax = self.general.n
        for i in range(len(self.STem)):
            if 0<n and n <= self.nmax:
                self.Snq1.append(self.STem[i])
            else:
                self.Snq1.append(0)
            
        self.Snq = np.array(self.Snq1).reshape(self.general.ntot , self.general.n)         
        self.rmn0=[]
        self.rmn1 = self.VCM.MatrixCreate(self.rmn0 , self.general.n , self.general.ntot)
        self.Rmnpq = Matrix(self.general.n*self.general.ntot , self.general.n*self.general.ntot)
        " this for loop changes M2Tem indices with function M2Tem by same order of (m,n,p,q) in m02"
        "Attention: m , n start from 1 but p,q start from 0"
        for i in range (1,self.general.n*self.general.ntot+1 ):
            for j in range (1,self.general.n*self.general.ntot+1 ):
                m = self.m2[i-1,j-1,0]
                n = self.m2[i-1,j-1,1]
                p = self.m2[i-1,j-1,2]-1
                q = self.m2[i-1,j-1,3]-1
                R = self.ab.smnpq(m,n,p,q)
                Matrix.setElement(self.Rmnpq,i,j,R)

        self.Rmnpq_array0 = []
        for i in range(self.general.n*self.general.ntot):
            for j in range(self.general.n*self.general.ntot):
                self.Rmnpq_array0.append(self.Rmnpq[i][j])
        
        self.Rmnpq_array = np.array(self.Rmnpq_array0).reshape(self.general.n*self.general.ntot ,self.general.n*self.general.ntot )
        self.f1 = np.array(np.dot(self.Rmnpq_array , self.STem)).reshape(self.general.ntot , self.general.n)
        self.fnq1 = np.array(self.f1).reshape( self.general.ntot , self.general.n)
        self.Tfar = []
        for q in range (self.general.ntot):
            self.Tfar.append(np.real(self.OmegaFar(q)))
        self.Tmatrix = []
        for q in range (self.general.ntot):
            self.Tmatrix.append(np.real(self.OmegaFar(q)+self.OmegaSQ(q)))
        self.TQ = []
        for q in range (self.general.ntot):
            self.TQ.append(np.real(self.OmegaGQ(q)))
        self.Tfieldtotal = self.TField(self.general.n , self.general.ntot) + np.conj(self.gamma)*self.z + self.t0
        self.Bettal0 = 2*self.TE.ThermalExpansion0 #or 2*self.TE.ThermalExpansion0*(1+self.general.Nu0)
        self.BettalQ = []
        for q in range(self.general.ntot):
            self.BettalQ.append(2*self.TE.ThermalExpansionQ[q])
        pass
    def KroneckerDelta(self,i,j):
        if i==j :
            x = 1
        else:
            x = 0
        return x
    def Fn(self,n,q):
        f = self.KroneckerDelta(n,1)*self.general.dd(q)*np.conj(self.gamma)/2
        return f

    def BCTem1(self,n,q):
        f = -(self.TE.ThermalConductivityQBar[q]-1)*(self.Fn(n,q)+np.conj(self.Fn(n,q))*self.general.v0(q)**(2*n))
        return f
    def M1Tem1(self,m,n,p,q):    
        f1 = (self.TE.ThermalConductivityQBar[q]-1)*self.ab.smnpq(m,n,p,q)
        f2 = (self.TE.ThermalConductivityQBar[q]+(self.general.v0(q)**(2*n) +self.general.v0(q)**(-2*n))/((self.general.v0(q)**(2*n) - self.general.v0(q)**(-2*n))))*self.KroneckerDelta(m,n)*self.KroneckerDelta(p,q)
        return f1+f2
    def M2Tem1(self,m,n,p,q):
        f1 = (self.TE.ThermalConductivityQBar[q]-1)*self.general.v0(q)**(2*n)*np.conj(self.ab.smnpq(m,n,p,q))
        f2 = (2/(self.general.v0(q)**(2*n) - self.general.v0(q)**(-2*n)))*self.KroneckerDelta(m,n)*self.KroneckerDelta(p,q)
        return f1+f2
    
    def fnq(self,n,q):

        if n>0 and n<=self.maxn:
            x = self.Fn(n,q) + self.fnq1[q][n-1]
        elif n<0 and n>=-self.maxn:
            x = self.Fn(-n,q) + self.fnq1[q][-n-1]
        elif n==0:
            x = self.Omegafar0[q]
        else:
            x=0
        return x
    def fnqr(self,n,q):
        if n>0 and n<=self.maxn:
            x = self.f1[q][n-1]
        elif n<0 and n>=-self.maxn:
            x = self.f1[q][-n-1]
        else:
            x=0
        return x
    def Gnq2(self,m,q):
        TxQ=self.TE.ThermalConductivityQBar[q] 
        n = abs(m)
        Vp = self.general.v0(q)**(4*n) -1
        f1 = self.general.v0(q)**(2*n) - self.general.v0(q)**(-2*n)
        if f1 == 0:
            f2 = self.fnq(n,q)
        else:
            f2 = (self.fnq(n ,q)*Vp - self.general.v0(q)**(2*n)*np.conj(self.Snq[q][n-1]) - self.Snq[q][n-1])/Vp
        return f2/TxQ

    def Gnq4(self,m,q):
        n = abs(m)
        if n == 0:
            f = self.fnq(n,q)
        else:
            f = 0.5*(2*self.fnq(n,q) + (self.general.v0(q)**(-n))*((np.conj(self.Snq[q][n-1])-self.Snq[q][n-1])/((self.general.v0(q)**n - self.general.v0(q)**(-n))) + (np.conj(self.Snq[q][n-1]) + self.Snq[q][n-1])/((self.general.v0(q)**(n) + self.general.v0(q)**(-n) ))))
        return f

    def OmegaFar(self,q) : 
        v = cmath.exp(self.Xi[q])
        sum = 0
        for i in range (1,self.general.n+1):
            x= self.fnq(i,q)*(v**i)*v**(-(i*i-1)/4)
            sum = sum + x
        sum1 = 0
        for j in range (-self.general.n-1 , 1):
            x= self.fnq(j,q)*(v**j)
            sum1 = sum1 + x
        return sum + sum1
    def OmegaSQ(self,q):
        v = cmath.exp(self.Xi[q])
        sum = 0
        for i in range (self.general.n):
            x= self.Snq[q][i]*(v**-i)
            sum = sum + x
        return sum

    def Gnq (self,m,q):
        x = complex(np.real(self.Gnq4(m,q)) , np.imag(self.Gnq2(m,q)))
        return x

    def OmegaGQ(self,q):
        v = cmath.exp(self.Xi[q])
        sum = 0
        for i in range (-self.general.n , self.general.n):
            x= self.Gnq(i,q)*(v**i)
            sum = sum + x
        return sum
    def TField(self , n , ntot):
        def sumq(q):
            sumq = 0
            for i in range (self.general.n):
                x= self.Snq[q][i]*cmath.exp(self.Xi[q])
                sumq = sumq + x
            return sumq
        sum0 = 0
        for q in range (ntot):
            x = sumq(q)*cmath.exp(complex( 0 , self.general.Theta[q]))
            sum0 = sum0 + x
        return sum0

    def bcTem1nq(self,n,q):
        f1 = self.Bettal0*self.general.dd(q)/n/2
        f2 = self.BettalQ[q]*self.general.dd(q)/n/2
        f3 = self.fnq(n+1,q) - self.fnq(n-1,q) - (self.Snq[q][n-2]-self.Snq[q][n]) - self.Snq[q][0]*(-self.general.v0(q))**n
        f4 = self.Gnq(n+1 , q) - self.Gnq(n-1,q)
        return f1*f3 + f2*f4

    def bcTem2nq(self , n,q):
        f0 = self.general.v0(q)**(2*n)
        f1 =self.Bettal0*self.general.dd(q)/n/2
        f2 = -self.BettalQ[q]*self.general.dd(q)/n/2
        f3 = (self.fnq(n-1,q) - self.fnq(n+1,q))*f0 - self.Snq[q][0]*(-self.general.v0(q))**n
        f4 = (self.Gnq(n-1,q) - self.Gnq(n+1,q))*f0
        return f1*f3+f2*f4
    




    
    


