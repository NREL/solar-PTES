# This is a python program to provide a simple model of Pumped Thermal Energy Storage (PTES)

# The PTES is a recuperated Joule-Brayton cycle with liquid storage tanks.
# See McTigue 2015, Laughlin 2019 for more details.

# Code developed by Josh McTigue, July 2020 for Jen King's NREL LDRD
# JoshuaDominic.McTigue@nrel.gov

# Module that specifies the input parameters
def inputs():

    inputs = {
        "eff": 0.97, # Heat exchanger effectiveness
        "eta": 0.90, # Polytropic efficiency of compressors and expanders
        "etaP":0.70, # Isentropic efficiency of pumps and fans
        "ploss": 0.01, # Fractional pressure loss in each heat exchanger (working fluid)
        "plossL":0.01, # Fractional pressure loss (liquid)
        "plossA":0.001, # Fractional pressure loss (air)
        "motor":0.98,   # Motor efficiency
        "gen":0.98,     # Generator efficiency
        "T0": 300,     # Ambient temperature, K
        "P0":1e5,      # Ambient pressure, Pa
        "P1": 10e5,    # Lowest pressure in the cycle, Pa
        "T1": 573,     # Charging compressor inlet temperature, K
        "T2": 843,     # Charging compressor outlet temperature, K
        "WF": "Nitrogen",     # Working fluid
        "HF": "NitrateSalt",  # Hot storage fluid
        "CF": "Glycol",       # Cold storage fluid
        "Pdis": 10e6,         # Power output, W
        "Tchg": 4,    # Charging duration, h
        "Tdis": 10,   # Discharging duration
        
    }

    return inputs

# This module calculates the properties around the charging cycle
def charge(inputs,cyc):

    # Unpack some classes first
    WF = cyc.WF
    HF = cyc.HF
    CF = cyc.CF
    A1 = cyc.A1

    # Unpack frequently used variables
    eff = cyc.eff
    eta = cyc.eta
    ploss = cyc.ploss
    plossL = cyc.plossL
    plossA = cyc.plossA 
    
    # Compressor inlet
    WF.Tc[0] = inputs['T1']
    WF.Tc[1] = inputs['T2']
    WF.Pc[0] = inputs['P1']
    WF.mdotC = 1 

    # Compressor outlet
    cyc.PRchgC = (WF.Tc[1] / WF.Tc[0]) ** (WF.gam * eta / (WF.gam-1))
    WF.Pc[1] = WF.Pc[0] * cyc.PRchgC

    # Have to go to recuperator first
    WF.Tc[6] = inputs['T0']
    WF.Tc[2] = (WF.Tc[0] + WF.Tc[6] * (eff-1))/eff
    WF.Tc[3] = WF.Tc[2] - WF.Tc[0] + WF.Tc[6]

    # Now try hot heat exchanger
    HF.Tc[0] = WF.Tc[0]
    HF.Tc[1] = HF.Tc[0] + eff * (WF.Tc[1] - HF.Tc[0])
    HF.mdotC = WF.mdotC*WF.cp*(WF.Tc[1]-WF.Tc[2]) / (HF.cp*(HF.Tc[1]-HF.Tc[0]))
    HF.Pc[0] = inputs['P0']
    HF.Pc[1] = HF.Pc[0] * (1 - plossL)

    # Sort out pressures
    WF.Pc[2] = WF.Pc[1] * (1 - ploss)
    WF.Pc[3] = WF.Pc[2] * (1 - ploss)
    
    # Hot heat exchanger
#    HF.Tc[0] = WF.Tc[0]
#    WF.Pc[2] = WF.Pc[1] * (1 - ploss)
#    WF.Tc[2] = WF.Tc[1] - 0.992*eff * (WF.Tc[1] - HF.Tc[0])
#    HF.Tc[1] = HF.Tc[0] + WF.Tc[1] - WF.Tc[2]

#    HF.mdotC = WF.mdotC*WF.cp*(WF.Tc[1]-WF.Tc[2]) / (HF.cp*(HF.Tc[1]-HF.Tc[0]))
#    HF.Pc[0] = inputs['P0']
#    HF.Pc[1] = HF.Pc[0] * (1 - plossL)

    # Recuperator
#    WF.Pc[3] = WF.Pc[2] * (1 - ploss)
#    WF.Tc[6] = (eff * WF.Tc[2] - WF.Tc[0]) / (eff - 1)
#    WF.Tc[3] = WF.Tc[2] - WF.Tc[0] + WF.Tc[6]

    # Heat rejection
    WF.Pc[4] = WF.Pc[3] * (1 - ploss)

    A1.Tc[0] = inputs['T0']
    WF.Tc[4] = WF.Tc[3] - eff * (WF.Tc[3] - A1.Tc[0])
    alp      = 2 # This is a ratio of mdot*cp
    A1.Tc[1] = A1.Tc[0] + (WF.Tc[3] - WF.Tc[4])/alp

    A1.mdotC = WF.mdotC*WF.cp*(WF.Tc[3]-WF.Tc[4])/(A1.cp*(A1.Tc[1]-A1.Tc[0]))
    A1.Pc[0] = inputs['P0']
    A1.Pc[1] = A1.Pc[0] * (1 - plossA)

    # Expansion
    cyc.PRchgE = cyc.PRchgC * (1-ploss)**(5)
    WF.Pc[5] = WF.Pc[4] / cyc.PRchgE

    WF.Tc[5] = WF.Tc[4] / (cyc.PRchgE ** ((WF.gam-1)*eta/WF.gam))

    # Cold heat exchanger - working fluid temperatures already found
    WF.Pc[6] = WF.Pc[5] * (1 - ploss)
    CF.Tc[0] = WF.Tc[5] + (WF.Tc[6]-WF.Tc[5])/eff
    CF.Tc[1] = CF.Tc[0] - WF.Tc[6] + WF.Tc[5]

    CF.mdotC = WF.mdotC*WF.cp*(WF.Tc[6]-WF.Tc[5])/(CF.cp*(CF.Tc[0]-CF.Tc[1]))
    CF.Pc[0] = inputs['P0']
    CF.Pc[1] = CF.Pc[0] * (1 - plossL)
    
    print('Design charging calculations complete.')

    print('\nWorking fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(WF.Tc)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,WF.Tc[i]-273,WF.Pc[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(WF.mdotC))

    print('\nHot storage fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(HF.Tc)):
        print('#{0:1d} {1:14.1f} {2:9.2f}'.format(i+1,HF.Tc[i]-273,HF.Pc[i]/1e5))
    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(HF.mdotC))

    print('\nCold storage fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(CF.Tc)):
        print('#{0:1d} {1:14.1f} {2:9.2f}'.format(i+1,CF.Tc[i]-273,CF.Pc[i]/1e5))
    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(CF.mdotC))

    print('\nAir stream 1 conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(A1.Tc)):
        print('#{0:1d} {1:14.1f} {2:9.2f}'.format(i+1,A1.Tc[i]-273,A1.Pc[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(A1.mdotC))   


    # Pack the classes back up
    cyc.WF = WF
    cyc.HF = HF
    cyc.CF = CF
    cyc.A1 = A1

    return cyc



# This module calculates the properties around the discharging cycle
def discharge(inputs,cyc):

    # Unpack some classes first
    WF = cyc.WF
    HF = cyc.HF
    CF = cyc.CF
    A1 = cyc.A1
    A2 = cyc.A2

    # Unpack frequently used variables
    eff = cyc.eff
    eta = cyc.eta
    ploss = cyc.ploss
    plossL = cyc.plossL
    plossA = cyc.plossA 


    # Index 0 is the compressor inlet
    WF.Pd[0] = inputs['P1']
    WF.mdotD = 1.0

    # ... but we start at the hot heat exchanger
    HF.Td[0] = HF.Tc[1]
    HF.Td[1] = HF.Tc[0] # Hot fluid must return to its original temperature

    # Mass flow is the same as during charge
    HF.mdotD = HF.mdotC
    mcpMIN   = min(HF.mdotD*HF.cp , WF.mdotD*WF.cp)
    WF.Td[3] = HF.Td[0] - HF.mdotD * HF.cp * (HF.Td[0] - HF.Td[1]) / (eff * mcpMIN)
    WF.Td[4] = WF.Td[3] + HF.mdotD * HF.cp * (HF.Td[0] - HF.Td[1]) / (WF.mdotD * WF.cp)
    
#    WF.Td[3] = HF.Td[0] - (HF.Td[0] - HF.Td[1]) / eff
#    WF.Td[4] = WF.Td[3] + HF.Td[0] - HF.Td[1]

#    HF.mdotD = WF.mdotD*WF.cp*(WF.Td[4]-WF.Td[3])/(HF.cp*(HF.Td[0]-HF.Td[1]))
    HF.Pd[0] = inputs['P0']
    HF.Pd[1] = HF.Pd[0] * (1 - plossL)

    # Recuperator
    WF.Td[2] = WF.Tc[4] # An assumption
    WF.Td[5] = WF.Td[2] + (WF.Td[3] - WF.Td[2])/eff
    WF.Td[6] = WF.Td[5] - WF.Td[3] + WF.Td[2]

    # Expander
    cyc.PRdisE = (WF.Td[4]/WF.Td[5]) ** (WF.gam / (eta * (WF.gam-1)))
    WF.Pd[5] = WF.Pd[0] / (1-ploss)**3
    WF.Pd[4] = WF.Pd[5] * cyc.PRdisE

    # Cold heat exchanger
    CF.Td[0] = CF.Tc[1]
    CF.Td[1] = CF.Tc[0] # Cold fluid returns to its original temperature
    WF.Td[7] = CF.Td[0] + (CF.Td[1] - CF.Td[0])/eff
    WF.Td[0] = WF.Td[7] - CF.Td[1] + CF.Td[0]

    CF.mdotD = WF.mdotD*WF.cp*(WF.Td[7]-WF.Td[0])/(CF.cp*(CF.Td[1]-CF.Td[0]))
    CF.Pd[0] = inputs['P0']
    CF.Pd[1] = CF.Pd[0] * (1 - plossL)
    
    # Compressor
    cyc.PRdisC = cyc.PRdisE / (1-ploss)**6
    WF.Td[1]   = WF.Td[0] * cyc.PRdisC ** ((WF.gam-1)/(WF.gam * eta))
    WF.Pd[1]   = WF.Pd[0] * cyc.PRdisC

    # Heat rejection 1
    A1.Td[0] = inputs['T0']
    alp = 2 # Ratio of mdot*cp
    A1.Td[1] = A1.Td[0] + eff * (WF.Td[1] - A1.Td[0])/alp

    A1.mdotD = WF.mdotD*WF.cp*(WF.Td[1]-WF.Td[2])/(A1.cp*(A1.Td[1]-A1.Td[0]))
    A1.Pd[0] = inputs['P0']
    A1.Pd[1] = A1.Pd[0] * (1 - plossA)    

    # Heat rejection 2
    A2.Td[0] = inputs['T0']
    A2.Td[1] = A2.Td[0] + eff * (WF.Td[6] - A1.Td[0])/alp
    
    A2.mdotD = WF.mdotD*WF.cp*(WF.Td[6]-WF.Td[7])/(A2.cp*(A2.Td[1]-A2.Td[0]))
    A2.Pd[0] = inputs['P0']
    A2.Pd[1] = A2.Pd[0] * (1 - plossA)

    # Calculate the remaining pressures
    WF.Pd[2] = WF.Pd[1] * (1-ploss)
    WF.Pd[3] = WF.Pd[2] * (1-ploss)
    WF.Pd[6] = WF.Pd[5] * (1-ploss)
    WF.Pd[7] = WF.Pd[6] * (1-ploss)
    

    print('\nDesign discharging calculations complete.')

    print('\nWorking fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(WF.Td)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,WF.Td[i]-273,WF.Pd[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(WF.mdotD))

    print('\nHot storage fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(HF.Td)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,HF.Td[i]-273,HF.Pd[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(HF.mdotD))
        
    print('\nCold storage fluid conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(CF.Td)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,CF.Td[i]-273,CF.Pd[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(CF.mdotD))
        
    print('\nAir stream 1 conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(A1.Td)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,A1.Td[i]-273,A1.Pd[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(A1.mdotD))
        
    print('\nAir stream 2 conditions:')
    print('State #     Temp, C    Pressure, bar')
    for i in range(len(A2.Td)):
        print('#{0:1d} {1:14.1f} {2:10.2f}'.format(i+1,A2.Td[i]-273,A2.Pd[i]/1e5))

    print('Relative mass flow rate (kg/s per kg/s working fluid):   {0:5.3f}'.format(A2.mdotD))    

        
    # Pack the classes back up
    cyc.WF = WF
    cyc.HF = HF
    cyc.CF = CF
    cyc.A1 = A1
    cyc.A2 = A2

    return cyc



# This class contains data for the PTES cycle
class Cycle:
    def __init__(self,inputs):
        self.eff = inputs['eff']
        self.eta = inputs['eta']
        self.etaP = inputs['etaP']
        self.ploss = inputs['ploss']
        self.plossL = inputs['plossL']
        self.plossA = inputs['plossA']
        self.motor  = inputs['motor']
        self.gen    = inputs['gen']
        
        self.Win   = 0 # Design work input (charge)
        self.Wout  = inputs['Pdis'] # Design work output (discharge)
        self.Ein   = 0 # Total electricity into the system, J
        self.Eout  = 0 # Total electricity out of the system, J
        self.Tdis  = inputs['Tdis'] * 3600
        self.Tchg  = inputs['Tchg'] * 3600

        self.PRchgC = 0 # Charging compression pressure ratio
        self.PRchgE = 0 # Charging expansion ratio
        self.PRdisC = 0 # Discharging compression pressure ratio
        self.PRdisE = 0 # Discharging expansion ratio
        self.win   = 0 # Specific design work input
        self.wout  = 0 # Specific design work output
        self.qin   = 0 # Heat into hot storage, charge
        self.qout  = 0 # Heat out of hot storage, discharge
        self.RTeff = 0 # Design round-trip efficiency
        self.HEeff = 0 # Heat engine efficiency
        self.COP   = 0 # Coefficient of performance
        
        self.mdotC = 0 # Design charging mass flow rate, kg/s
        self.mdotD = 0 # Design discharging mass flow rate, kg/s

        # Set up fluid classes
        self.WF = Fluid(inputs['WF'],'WF')  # Working fluid
        self.HF = Fluid(inputs['HF'],'HF') # Hot fluid
        self.CF = Fluid(inputs['CF'],'CF') # Cold fluid
        self.A1 = Fluid('Air','Air')
        self.A2 = Fluid('Air','Air')
        self.A3 = Fluid('Air','Air')
        print('Set up fluid classes')


    # Calculate the specific work input and output etc.    
    def performance (self):

        # Unpack some classes first
        WF = self.WF
        HF = self.HF
        CF = self.CF
        A1 = self.A1
        A2 = self.A2

        Tc = WF.Tc
        Td = WF.Td

        # Calculate the specific works
        # Charge
        self.win = WF.cp * ((Tc[1]-Tc[0]) - (Tc[4]-Tc[5]))
        HFpara   = (HF.mdotC / HF.rho) * (HF.Pc[0] - HF.Pc[1]) / self.etaP
        CFpara   = (CF.mdotC / CF.rho) * (CF.Pc[0] - CF.Pc[1]) / self.etaP

        pr = 1 / (1-self.plossA)
        A1para = A1.mdotC*A1.cp*A1.Tc[0]*(pr**((A1.gam-1)/A1.gam)-1)/self.etaP
#        A1para   = (A1.mdotC / A1.rho) * (A1.Pc[0] - A1.Pc[1]) / self.etaP

        self.win = self.win + HFpara + CFpara + A1para
        self.win = self.win / self.motor

        # Discharge
        self.wout = WF.cp * ((Td[4]-Td[5]) - (Td[1]-Td[0]))
        HFpara   = (HF.mdotD / HF.rho) * (HF.Pd[0] - HF.Pd[1]) / self.etaP
        CFpara   = (CF.mdotD / CF.rho) * (CF.Pd[0] - CF.Pd[1]) / self.etaP

        pr = 1 / (1-self.plossA)
        A1para = A1.mdotD*A1.cp*A1.Td[0]*(pr**((A1.gam-1)/A1.gam)-1)/self.etaP
        A2para = A1.mdotD*A1.cp*A1.Td[0]*(pr**((A1.gam-1)/A1.gam)-1)/self.etaP

        self.wout = self.wout - HFpara - CFpara - A1para - A2para
        self.wout = self.wout * self.gen
        
        self.qin  = WF.cp * (Tc[1]-Tc[2])
        self.qout = WF.cp * (Td[4]-Td[3])
        
        self.RTeff = self.wout / self.win
        self.COP   = self.qin / self.win
        self.HEeff = self.wout / self.qout

        self.Eout  = self.Wout * self.Tdis
        self.Ein   = self.Eout / self.RTeff
        self.Win   = self.Ein / self.Tchg

        self.mdotD = self.Wout / self.wout
        self.mdotC = self.Win / self.Tchg

        print('\nDesign performance calculations complete:')

        print('Design power input:  {0:10.1f} MW'.format(self.Win/1e6))
        print('Design power output: {0:10.1f} MW\n'.format(self.Wout/1e6))

        print('Design energy input:  {0:10.1f} MWh'.format(self.Ein/3600e6))
        print('Design energy output: {0:10.1f} MWh\n'.format(self.Eout/3600e6))

        print('Round-trip efficiency:      {0:10.2f} %'.format(self.RTeff * 100))
        print('Coefficient of Performance: {0:10.2f} '.format(self.COP))
        print('Heat engine efficiency:     {0:10.2f} %\n'.format(self.HEeff * 100))

        return self
        
        
    def off_design(self):
        print('Off-design behaviour is not yet implemented')

    def state_of_charge(self):
        print('State-of-charge calculations not yet implemented')

        
            
    

# This class contains fluid properties
class Fluid:
    def __init__(self,name,type):
        self.name = name
        self.type = type

        self.mdotC = 0
        self.mdotD = 0
        if self.type == "WF":
            if self.name == "Nitrogen":
                self.cp = 1041.3
                self.cv = 743.2
                self.rho = 1.123
                self.gam = self.cp / self.cv
            elif self.name == "Argon":
                self.cp = 521.5
                self.cv = 312.4
                self.rho = 1.603
                self.gam = self.cp / self.cv
            elif self.name == "Hydrogen":
                self.cp = 14312.8
                self.cv = 10186.4
                self.rho = 0.0807
                self.gam = self.cp / self.cv
            elif self.name == "Helium":
                self.cp = 5193.2
                self.cv = 3116.1
                self.rho = 0.1603
                self.gam = self.cp / self.cv
            else:
                err = self.name+" is not a valid fluid for the working fluid"
                raise NameError(err)

            self.Tc = [-1 for x in range(7)]
            self.Pc = [-1 for x in range(7)]

            self.Td = [-1 for x in range(8)]
            self.Pd = [-1 for x in range(8)]

            

        elif self.type == "HF":
            if self.name == "NitrateSalt":
                self.cp = 2000
                self.rho = 800
                self.Tmax = 570 + 273
                self.Tmin = 250 + 273
            elif self.name == "ChlorideSalt":
                self.cp = 2000
                self.rho = 800
                self.Tmax = 750 + 273
                self.Tmin = 400 + 273
            else:
                err = self.name+" is not a valid fluid for the hot storage fluid"
                raise NameError(err)

            self.Tc = [-1 for x in range(2)]
            self.Pc = [-1 for x in range(2)]

            self.Td = [-1 for x in range(2)]
            self.Pd = [-1 for x in range(2)]


        elif self.type == "CF":
            if self.name == "Glycol":
                self.cp = 2000
                self.rho = 800
                self.Tmax = 200 + 273
                self.Tmin = -50 + 273
            elif self.name == "Methanol":
                self.cp = 2000
                self.rho = 800
                self.Tmax = 200 + 273
                self.Tmin = -100 + 273
            else:
                err = self.name+" is not a valid fluid for the cold storage fluid"
                raise NameError(err)

            self.Tc = [-1 for x in range(2)]
            self.Pc = [-1 for x in range(2)]

            self.Td = [-1 for x in range(2)]
            self.Pd = [-1 for x in range(2)]

        elif self.type == "Air":
            self.cp = 1000
            self.cv = 718
            self.gam = self.cp / self.cv
            self.rho = 1.225
            self.Tc = [-1 for x in range(2)]
            self.Pc = [-1 for x in range(2)]

            self.Td = [-1 for x in range(2)]
            self.Pd = [-1 for x in range(2)]
