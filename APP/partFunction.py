##############################################################

# -*- coding: utf-8 -*-

__author__ = "Thiago Lopes, Mateus Barbosa and Daniel Machado"
__credits__ = ["LMSC"]
__version__ = "1.0.0"
__maintainer__ = "Thiago Lopes"
__email__ = "lopes.th.o@gmail.com"
__status__ = "Development"

# Description: 

# Built-in/Generic Imports

# Libs

# Own modules
from numpy import pi, power, exp, arange
import matplotlib.pyplot as plt

# Constant Values
N_a = 6.02214076E23
c = 299792458
h = 6.62607015E-34

def getIdealGasConstant(unit = "j/K.mol"):
    if unit == "j/K.mol":
        return 8.31446261815324
    elif unit ==  "m3.Pa/K.mol":
        return 8.31446261815324
    elif unit == "L.bar/K.mol":
        return 8.31446261815324E-2

class ReducedMass(object):
    
    def __init__(self, m1, m2):
        self.__mass = (m1 * m2) / (m1 + m2)
        
    def getReducedMass(self):
        return self.__mass

class Qtrans(object):
    
    def __init__(self, temperatures, mass):
        self.__temperatures = temperatures
        self.__mass = mass
        self.__R = getIdealGasConstant("L.bar/K.mol")
        self.__q = []
        
    def values(self, pressure):
        if len(self.__q) == 0:
            self.__calc(pressure)
        return self.__q
        
    def __calc(self, p):
        for T in self.__temperatures:
            part1 = power(((2*pi*self.__mass*self.__R * T)/(power(N_a, 2)*power(h, 2))), 1.5)
            part2 = (self.__R * T)/(N_a * p)
            self.__q.append(part1 * part2)

class Qelet(object):

    def __init__(self, temperatures, multiplicity, De):
        self.__multiplicity = multiplicity
        self.__De = De
        self.__temperatures = temperatures
        self.__R = getIdealGasConstant("j/K.mol")
        self.__q = []
    
    def values(self):
        if len(self.__q) == 0:
            self.__calc()
        return self.__q
    
    def __calc(self):
        for T in self.__temperatures:
            self.__q.append(self.__multiplicity * exp((self.__De)/(self.__R * T)))    

class Qrot(object):
    
    def __init__(self, temperatures, reducedMass, Re):
        self.__temperatures = temperatures
        self.__Be = (h * N_a) / (8 * power(pi, 2) * c * reducedMass * power(Re, 2)) * 1.196266E3
        self.__R = getIdealGasConstant("j/K.mol")
        self.__q = []
    
    def values(self):
        if len(self.__q) == 0:
            self.__calc()
        return self.__q
        
    def __calc(self):
       for T in self.__temperatures:
           value = (self.__R * T / self.__Be) + (1/3) + (1/15)*(self.__Be / (self.__R *T)) + (4/315)*power((self.__Be / (self.__R *T)), 2) + (1/315)*power((self.__Be / (self.__R *T)), 3)
           self.__q.append(value)
        
class Qvib(object):
    
    def __init__(self, vMax, w_e, w_exe, w_eye, temperatures):
        self.__temperatures = temperatures
        self.__vMax = vMax
        self.__we = w_e
        self.__wexe = w_exe
        self.__weye = w_eye
        self.__R = getIdealGasConstant("j/K.mol")
        self.__q = []
    
    def values(self):
        if len(self.__q) == 0:
            self.__calc()
        return self.__q
    
    def __calc(self):
        for T in self.__temperatures:
            value = 0
            for v in arange(self.__vMax+1):
                value += exp(-self.__E_v(v)/(self.__R * T))
            self.__q.append(value)
    
    def __E_v(self, v):
        return self.__we * (v + 0.5) - self.__wexe * power((v + 0.5), 2) + self.__weye * power((v + 0.5), 3)


class PartitionFunction(object):
    
    def __init__(self, rot_N_vib = False):
        self.__rot_N_vib = rot_N_vib
        self.__q_tot = []
    
    def values(self, temperatures, p, multiplicity, De, m1, m2 = None, Re = None, vMax = None, w_e = None, w_exe = None, w_eye = None):

        if m2 != None:
            mi = ReducedMass(m1, m2).getReducedMass()
            m_t = m1 + m2
        else:
            mi = m1
            m_t = m1
        q_trans = Qtrans(temperatures, m_t).values(p)
        q_ele = Qelet(temperatures, multiplicity, De).values()
        if self.__rot_N_vib:
            q_rot = Qrot(temperatures, mi, Re).values()
            q_vib = Qvib(vMax, w_e, w_exe, w_eye, temperatures).values()
            
        for i in arange(len(temperatures)):
            value = q_trans[i] * q_ele[i]
            if self.__rot_N_vib:
                value = value * q_rot[i]
                value = value * q_vib[i]
            self.__q_tot.append(value)
        return(self.__q_tot)
    
    
if __name__ == "__main__":
    x = arange(90, 200, 1)
    
    '''
    y = PartitionFunction().values(x, 1E-6, 1, 0, 32.0419E-3)  # Função de Partição Total sem rotação e vibração   
    y = Qrot(x, ReducedMass(32.0419E-3, 39.948E-3).getReducedMass(), 2.86E-10).values()
    y = Qvib(vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3, temperatures = x).values()
    y = Qelet(temperatures = x, multiplicity = 1, De = 1.3252E3).values()
    y = Qtrans(x,32.0419E-3 + 39.948E-3).values(1E-6)
    '''
    y = PartitionFunction(True).values(temperatures = x, p =  1E-6, multiplicity = 1, De = 1.3252E3, m1 = 32.0419E-3, m2 = 39.948E-3, Re = 2.86E-10, vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3)
    print(y[0])
    p = plt.plot(x, y)
    plt.show()
    




