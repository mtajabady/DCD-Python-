from readInputFiles import RWfile
class ThermalExpansion:
    def __init__(self,TS):
        ThermalExpansion = RWfile.readFile("D:/Dr_Ebrahimi/DCD_Package-main/output/ThermalExpansion-2holes")

        self.ThermalExpansion0 = float(ThermalExpansion[0])
        self.ThermalExpansiontemp = RWfile.list_input_float(ThermalExpansion , 1)
        self.ThermalConductivity0 = float(ThermalExpansion[2])
        self.ThermalConductivitytemp =RWfile.list_input_float(ThermalExpansion , 3)
        self.ThermalExpansionQ = []
        self.ThermalConductivityQ =[]
        self.ThermalConductivityQBar = []
        if TS == 1 :

            self.ThermalExpansionQ= self.ThermalExpansiontemp
            self.ThermalConductivityQ= self.ThermalConductivitytemp
            if self.ThermalConductivity0 == 0:
                self.ThermalConductivity0 = 0.0001
            else:
                for x in self.ThermalConductivitytemp:
                    self.ThermalConductivityQBar.append(x/self.ThermalConductivity0)

        pass