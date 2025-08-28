import cmath
import math 
from decimal import *
getcontext().prec += 5  # Extra digits for intermediate steps.

D = Decimal

import numpy as np
from scipy.integrate import quad
from generalCoefs import inputCoefs

def makeDecimal (x) :
        
    if isinstance(x, Decimal) :
        return (x)
    if isinstance(x, int) or isinstance(x, float) :
        return (Decimal(str(x)))
    if isinstance(x, complex) :
        return ( Decimal(str(x.real)), Decimal(str(x.imag)) )
    if isinstance(x, str) :
        status = 0
        try : x1 = Decimal(x)
        except : status = -1
        if not status : return x1

        status = 0
        try : x1 = eval(x)
        except : status = -1
        if status :
            print ('makeDecimal: bad status after eval.')
            exit(79)

        if isinstance(x1, str) :
            print ('makeDecimal: eval must not produce another string.')
            exit(78)
        return makeDecimal(x1)
    print ('makeDecimal: input must be Decimal, int, float, complex or str.')
    exit(77)

def largfactorial(n):       #gives natural logarithm of n!    
    return  n*np.log(n) - n + 1
    
class EtaMu():
    
    def __init__(self,general:inputCoefs) :
        self.general = general
        pass
    def complex_quadrature(self,func, a, b, **kwargs):           # Integrates on complex function "func" = lambda x : function(x) from a to b 
        def real_func(x):
            return np.real(func(x))
        def imag_func(x):
            return np.imag(func(x))
        real_integral = quad(real_func, a, b, **kwargs)
        imag_integral = quad(imag_func, a, b, **kwargs)
        return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

    def cpower(self, complexNumbwer , power):
    
        b = cmath.polar(complexNumbwer)
        x = math.pow(b[0] , power) * math.cos(power*b[1]) 
        y = math.pow(b[0] , power) * math.sin(power*b[1]) 

        return complex(x,y)

    def zpq(self , p,q):
        return self.general.zcenter[q] - self.general.zcenter[p]

    def vpq(self ,p,q):
        f = (self.zpq(p,q)/self.general.dpq(p,q)) + cmath.sqrt((self.zpq(p,q)/self.general.dpq(p,q)) + 1)*(cmath.sqrt((self.zpq(p,q)/self.general.dpq(p,q))-1))
        return f

    #vpqprime = (1/(self.general.dpq(p,q)))*(1 +  ((zpq(p, q)/self.general.dpq(p, q))/(cmath.sqrt((zpq(p,q)/self.general.dpq(p,q)) + 1)*(cmath.sqrt((zpq(p,q)/self.general.dpq(p,q)) - 1)))))
    def divPart1 (self,k,l,m,n):             #   Denomerator in M_nml
        x = math.factorial(k)*math.factorial(l-k)*math.factorial(k+n)*math.factorial(abs(m) +l-k)
        return x

    

    def M_nml(self,n,m,p,q,l):                    #   calculation of M_nml for k = 0 to k = l
        
        f2 = self.general.dd(p)/self.general.dd(q)

        sum01 = 0
        for k in range(l+1):
           # x = f2**(2*k) / self.divPart1
            x=(f2**(k))/self.divPart1(k,l,m,n)
            y = x*(f2**k)
            sum01 = sum01 + y
        return (sum01)

    def sum_lj(self,n,m,p,q,j):                      #   Sum from l to J in EtaClosePQ & MuClosePQ
        sumlj = 0
        d2 = self.general.dd(q)/self.general.dpq(p, q)
        for l in range(j+1):
            x=((-1)** (j - l)/math.gamma(j - l + 1))*(d2**(2* l + m) )*(self.M_nml(n,m,p,q,l))*(math.gamma(n + m + l + j ))
            sumlj = sumlj + x
        
        return sumlj

    " to avoid overflow error we used Ln in calculations"
    def sumFactor(self,n,m,l,p,q,xmu):       #   final Sum fator to use in total_sum
        fa = (n+abs(m) + 2*l - 1 + xmu)
        if fa < 20:
            a = np.log(math.factorial(fa))
        else :
            a = largfactorial(fa)
        b = ( n + abs(m)+ 2*l + xmu)*np.log(2*self.zpq(p,q))
        x = cmath.exp(a - b)
        return x

    def total_sum(self , n,m,p,q,degree,xmu):          #   calculation of second sigma for l = 0 to l =  degree
        sum02 = 0
        for l in range (degree+1):
            x= (self.general.dd(q)**(2*l + abs(m))) * self.M_nml(n,m,p,q,l) * self.sumFactor(n,m,l,p,q,xmu)
            sum02 = sum02 + x
        return sum02 

    def f1_calculations(self,n,m,p,q,xmu):           #   Xmu is 0 in Etafar and 1 in Mufar                        
        x = (1-3*xmu)*(n*(self.general.dd(p)** n)*((-1)** m) )*(self.total_sum(n,m,p,q,90,xmu)) #..........final formula for Etafar and Mufar
        return x 

    def sum_CLs(self,n,m,p,q,xmc):                   #   Total Sum in EtaClosePQ & MuClosePQ

        Sum_Cls = 0
        for j in range(19):
            if xmc == 0:
                a = 1
            else:
                a= ((n + m + 2*j))

            x= a*(self.vpq(p, q)**-(n + m + 2*j + xmc)) * self.sum_lj(n,m,p,q,j)
            Sum_Cls = Sum_Cls + x
        return Sum_Cls

    def f1_Cls(self,n,m,p,q,xmc):                    #   Final Function  for EtaClosePQ & MuClosePQ calculation : xmc is 0 in Eta and 1 in Mu 
        if xmc == 0:
            a=1
        else: 
            vpqprime = (1/(self.general.dpq(p,q)))*(1 +  ((self.zpq(p, q)/self.general.dpq(p, q))/(cmath.sqrt((self.zpq(p,q)/self.general.dpq(p,q)) + 1)*(cmath.sqrt((self.zpq(p,q)/self.general.dpq(p,q)) - 1)))))
            a = -vpqprime   
        d1 = self.general.dd(q)/self.general.dpq(p, q)
        x = (a*(n *(d1**n)*((-1)**m)))*self.sum_CLs(n,m,p,q,xmc)
        return x
   
    #(((EtafarPQ)))...........
    def EtafarPQ(self,n,m,p,q):
        if p==q:
            x = 0
        else:
            x = self.f1_calculations(n,m,p,q,0)
        return x

    #(((MufarPQ)))
    def MufarPQ(self ,n,m,p,q):
        if p==q:
            x = 0
        else:
            x = self.f1_calculations(n,m,p,q,1)
        return x

    #(((EtaClosePQ)))

    def EtaClosePQ(self,n,m,p,q):
        if p==q:
            x = 0
        else:
            x = self.f1_Cls(n,m,p,q,0)
        return x


    #(((MuClosePQ)))
    def MuClosePQ(self,n,m,p,q):
        if p==q:
            x = 0
        else:
            x = self.f1_Cls(n,m,p,q,1)
        return x

    def EtaIntgPQ (self,n,m,p,q):   
        if p==q:
            x=0
        else:
            a = self.complex_quadrature(lambda x: cmath.exp(-n*(cmath.acosh((self.general.dd(q)*math.cos(x)+self.zpq(p,q))/self.general.dd(p))))*math.cos(m*x) , 0, math.pi)       
            b = [np.real(a[0])/math.pi , np.imag(a[0])/math.pi]     
            x = complex(b[0] , b[1])
        return x
   
    def vn(self,z,n,p):    # when lim (vn -> 1) we have division by 0 in MuIntPQ so we have used "limit" Parameter in MuIntPQ 
        x = cmath.exp(-n*(cmath.acosh(z/self.general.dd(p))))
        return x
    #..........MuIntPQ
    
    def MuIntgPQ (self,n,m,p,q):
        
        if p==q:
            x=0
        else:
            if self.vn(self.general.dd(q)*math.cos(0)+self.zpq(p,q),n,p) == 1:
                limit = 0.00000589
            else:
                limit = 0 
            a =  self.complex_quadrature(lambda x: cmath.exp(-n*(cmath.acosh((self.general.dd(q)*math.cos(x)+self.zpq(p,q))/self.general.dd(p))))*math.cos(m*x)/(cmath.exp(-1*(cmath.acosh((self.general.dd(q)*math.cos(x)+self.zpq(p,q))/self.general.dd(p))))- (1/(cmath.exp(-1*(cmath.acosh((self.general.dd(q)*math.cos(x)+self.zpq(p,q))/self.general.dd(p))))))) , limit, math.pi)
            general = inputCoefs()
            x = 2*n*a[0]/(general.dd(p)*math.pi)

        return x
