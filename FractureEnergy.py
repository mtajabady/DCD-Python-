from readInputFiles import RWfile
class FractureEnergy:
    def __init__(self):
                        
        Fracture_Energy = RWfile.readFile("D:/Dr_Ebrahimi/DCD_Package-main/output/FractureEnergy")
        self.G0 = float(Fracture_Energy[0])
        self.GQ = RWfile.list_input_float(Fracture_Energy , 1)
        self.kQ=  RWfile.list_input_float(Fracture_Energy , 2)
        self.GGB= RWfile.list_input_float(Fracture_Energy , 3)
        self.vertexlist = []
        vt= RWfile.list_input_float(Fracture_Energy , 4)
        for i  in range (0,len(vt),2):
            self.vertexlist[len(self.vertexlist):] = [complex(vt[i],vt[i+1])]
        self.BE = int(Fracture_Energy[5])
        self.GB = int(Fracture_Energy[6])
        self.TS = int(Fracture_Energy[7])
        self.ltemp = float(Fracture_Energy[8])
        self.lnew = float(Fracture_Energy[9])
        self.delta = float(Fracture_Energy[10])
        self.l2ratio = float(Fracture_Energy[11])
        self.ttest = int(Fracture_Energy[12])
        self.inc = int(Fracture_Energy[13])
        self.naccuracy = int(Fracture_Energy[14])
    def printContent(self):
        print(self.GB)