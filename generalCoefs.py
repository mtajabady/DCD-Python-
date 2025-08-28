from readInputFiles import RWfile
import math , cmath
class inputCoefs:
    def __init__(self):
        '''    
        Geometrical properties of inclusions

        Parameters
        ----------
        n : int
            The total number of inclusions in the matrix.
        num1 : int
            Crack propagation criterion. Number 1, is the to search on the surface of inclusion for max energy, number 2 uses SIFs.
        num2: int
            It determines which process has to be used to calculate SIF; 1 is the formula from MM, 2 is the average SIF over short distance.
        Nu0: Float :
            Nu0 matrix poisson ratio.
        Float Mu0::usage ="Mu0 matrix shear modulus";

        (*Far-Field*)
        Float s11::usage ="s11 Farfield stress field";
        Float s12::usage ="s12 Farfield stress field";
        Float s22::usage ="s22 Farfield stress field";

        (*Properties of Inclusions*)
        Float Nu::usage ="Poisson ratio for the inclusions (ntot terms, ntot = number of inclusions).";
        Float MuQ::usage ="Shear modulus for the inclusions (ntot terms. Read information about the inclusion from the input file.";
        Float l1::usage ="l1=major axis of elliptical inclusions (ntot terms).";
        Float l2::usage ="l2=Minor axis of elliptical inclusions (ntot terms).";
        Float Theta::usage ="angle of inclination for each inclusion [-90,90]";
        Float Thetaprim::usage ="angle of inclination for each inclusion [0,360]";
        Complex zcenter::usage=" (x+Iy): center point of new inclusions(cracks) in the global coordinate system"; 
        Complex zcenternewlist1=" center point of new inclusions(cracks) in the global coordinate system"; 

            
            
            '''

        
        inputs = RWfile.readFile("D:/Dr_Ebrahimi/DCD_Package-main/output/output")
        self.n = int(inputs[0])
        self.nmax = 10    
        self.num1 = int(inputs[1])
        self.num2 = int(inputs[2])
        self.ntot = int(inputs[3])
        assert (self.ntot % 2) == 1, "ntot should be odd!"
        self.Nu0 = float(inputs[4])
        self.Mu0 = float(int(inputs[5]))
            
        #(*Far-Field*)
        self.s11 = float(inputs[6])
        self.s12 = float(inputs[7])
        self.s22 = float(inputs[8])

        #(*Properties of Inclusions*)
        
        self.Nu = RWfile.list_input_float(inputs, 9)
        self.MuQ = RWfile.list_input_float(inputs, 10)
        self.l1= RWfile.list_input_float(inputs, 11)
        self.l2 = RWfile.list_input_float(inputs, 12)
        self.Theta = RWfile.list_input_float(inputs, 13)
        for j in range(len(self.Theta)):
            self.Theta[j] = self.Theta[j]*math.pi/180.0

        self.Thetaprim = RWfile.list_input_float(inputs, 14)

        for j in range(len(self.Thetaprim)):
            self.Thetaprim[j] = self.Thetaprim[j]*math.pi/180.0

        self.zcenter = RWfile.list_input_complex(inputs, 15)
        self.zcenternewlist = RWfile.list_input_complex(inputs, 16)
        
        self.Chiq =[]
        for j in range(len(self.Nu)):
            chiq_0 = (3 - self.Nu[j])/(1+ self.Nu[j])
            chiq_1 = 3-4*self.Nu[j]
            self.Chiq.append([chiq_0 , chiq_1])

        self.Chi0 = 3 - 4*self.Nu0; 
            #(3 - Nu0)/(1+ Nu0);(*3 - 4*Nu0 This is for Plane strain*)

        self.MuBar = []
        for q in range(len(self.MuQ)):
            MuBar_p = (self.MuQ[q]/ self.Mu0)
            self.MuBar.append(MuBar_p)

        self.e = []
        for i in range(len(self.l2)):
            e_p = (self.l2[i]/ self.l1[i])
            self.e.append(e_p)

        self.E0 = 2.*self.Mu0*(1+self.Nu0)
        self.EQ = []
        for i in range(len(self.MuQ)):
            EQx = 2*self.MuQ[i]*(1+self.Nu[i])
            self.EQ.append(EQx)

                
        self.EndPoint = []
        for i in range(len(self.l1)):
            endpoint_i = self.zcenter[i]+self.l1[i]*cmath.exp(complex(0 , self.Thetaprim[i]))
            self.EndPoint.append(endpoint_i) 

        self.StartPoint=[]
        for i in range(len(self.l1)):
            startpoint_i = self.zcenter[i]-self.l1[i]*cmath.exp(complex(0 , self.Thetaprim[i]))
            self.StartPoint.append(startpoint_i)
    
    def Dq(self,p):
        x= math.sqrt(math.pow(self.l1[p],2) - math.pow(self.l2[p],2))
        return x

    def dd(self,p):
        f1 = self.Dq(p)
        f2 = cmath.exp(complex(0,self.Theta[p]))*f1
        return  f2  

    def dpq(self,p,q):
        x = self.dd(p)+self.dd(q)
        return x

    def Zeta0(self,p):
        if self.e[p] == 1:
            f = 0
        else:
            f = 0.5*math.log((1 + self.e[p])/(1 - self.e[p]))
        return f

    def v0(self,p) :
        f = math.exp(self.Zeta0(p))
        return f

    def fTheta(self,p ,q ):
        f= self.Theta[q] - self.Theta[p]
        return f

        


