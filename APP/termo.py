from numpy import array, diff, log, power, arange

from APP.partFunction import *

from matplotlib.pyplot import plot, show

class InternalEnergy(object):
    
    def __init__(self, temperatures, pressure, rot_N_vib = False):
        self.__R = getIdealGasConstant("j/K.mol")
        self.__rot_N_vib = rot_N_vib
        self.__t = temperatures
        self.__energy = []
        self.__p = pressure

    def setParameters(self, multiplicity, De, m1, m2 = None, Re = None, vMax = None, w_e = None, w_exe = None, w_eye = None):
        self.__multiplicity = multiplicity
        self.__De = De
        self.__m1 = m1
        self.__m2 = m2
        self.__Re = Re
        self.__vMax = vMax
        self.__w_e = w_e
        self.__w_exe = w_exe
        self.__w_eye = w_eye
        
    def values(self):
        if self.__multiplicity == None:
            print("Use the method InternalEnergy.setParameters(multiplicity, De, m1, m2(opt), Re(opt), vMax(opt), w_e(opt), w_exe(opt), w_eye(opt)) to indicate full partition function parameters")
            exit(0)
        self.__calc()
        return self.__energy
    
    def __calc(self):
        y = PartitionFunction(self.__rot_N_vib).values(self.__t, self.__p, self.__multiplicity, self.__De, self.__m1, self.__m2, self.__Re, self.__vMax, self.__w_e, self.__w_exe, self.__w_eye)
        diff_y = diff(log(y))
        temperatures = self.__t[1:]
        for i in arange(len(temperatures)):
            self.__energy.append((self.__R * power(temperatures[i], 2)) * diff_y[i])

class Enthalpy(InternalEnergy):
    
    def __init__(self, temperatures, pressure, rot_N_vib=False):
        super().__init__(temperatures, pressure, rot_N_vib=rot_N_vib)
        self.__H = []
        
    def values(self):
        if self._InternalEnergy__multiplicity == None:
            print("Use the method InternalEnergy.setParameters(multiplicity, De, m1, m2(opt), Re(opt), vMax(opt), w_e(opt), w_exe(opt), w_eye(opt)) to indicate full partition function parameters")
            exit(0)
        u = super().values()
        temperatures = self._InternalEnergy__t[1:]
        for i in arange(len(temperatures)):
            self.__H.append(u[i] + (self._InternalEnergy__R * temperatures[i]))
        return self.__H
            
class Entropy(InternalEnergy):
    
    def __init__(self, temperatures, pressure, rot_N_vib=False):
        super().__init__(temperatures, pressure, rot_N_vib=rot_N_vib)
        self.__energy = []
    
    def values(self):
        if self._InternalEnergy__multiplicity == None:
            print("Use the method InternalEnergy.setParameters(multiplicity, De, m1, m2(opt), Re(opt), vMax(opt), w_e(opt), w_exe(opt), w_eye(opt)) to indicate full partition function parameters")
            exit(0)
        self.__calc()
        return self.__energy        
    
    def __calc(self):
        y = PartitionFunction(self._InternalEnergy__rot_N_vib).values(self._InternalEnergy__t, self._InternalEnergy__p, self._InternalEnergy__multiplicity, self._InternalEnergy__De, self._InternalEnergy__m1, self._InternalEnergy__m2, self._InternalEnergy__Re, self._InternalEnergy__vMax, self._InternalEnergy__w_e, self._InternalEnergy__w_exe, self._InternalEnergy__w_eye)
        diff_y = diff(log(y))
        ln_y = log(y)[1:]
        temperatures = self._InternalEnergy__t[1:]
        for i in arange(len(temperatures)):
            self.__energy.append(self._InternalEnergy__R * temperatures[i] * diff_y[i] + (self._InternalEnergy__R * ln_y[i]))

class GibbsFreeEnergy(object):
    
    def __init__(self, temperatures, pressure, rot_N_vib=False):
        self.__R = getIdealGasConstant("j/K.mol")
        self.__rot_N_vib = rot_N_vib
        self.__t = temperatures
        self.__energy = []
        self.__p = pressure

    def setParameters(self, multiplicity, De, m1, m2 = None, Re = None, vMax = None, w_e = None, w_exe = None, w_eye = None):
        self.__multiplicity = multiplicity
        self.__De = De
        self.__m1 = m1
        self.__m2 = m2
        self.__Re = Re
        self.__vMax = vMax
        self.__w_e = w_e
        self.__w_exe = w_exe
        self.__w_eye = w_eye
        
    def values(self):
        if self.__multiplicity == None:
            print("Use the method InternalEnergy.setParameters(multiplicity, De, m1, m2(opt), Re(opt), vMax(opt), w_e(opt), w_exe(opt), w_eye(opt)) to indicate full partition function parameters")
            exit(0)
        self.__calc()
        return self.__energy
    
    def __calc(self):
        entropy = Entropy(self.__t, self.__p, self.__rot_N_vib)
        entropy.setParameters(self.__multiplicity, self.__De, self.__m1, self.__m2, self.__Re, self.__vMax, self.__w_e, self.__w_exe, self.__w_eye)
        entropy = entropy.values()
        enthalpy = Enthalpy(self.__t, self.__p, self.__rot_N_vib)
        enthalpy.setParameters(self.__multiplicity, self.__De, self.__m1, self.__m2, self.__Re, self.__vMax, self.__w_e, self.__w_exe, self.__w_eye)
        enthalpy = enthalpy.values()
        temperatures = self.__t[1:]
        for i in arange(len(temperatures)):
            self.__energy.append(enthalpy[i] - (temperatures[i]*entropy[i]))

if __name__ == "__main__":
    x = arange(90, 200, 1)
    '''
    teste = InternalEnergy(x, 1E-6, True)
    teste.setParameters(multiplicity = 1, De = 1.3252E3, m1 = 32.0419E-3, m2 = 39.948E-3, Re = 2.86E-10, vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3)
    y = teste.values()
    y = [value/1000 for value in y]
    print(y[0], y[-1])
    plot(x[1:], y)
    show()
    teste = Enthalpy(x, 1E-6, True)
    teste.setParameters(multiplicity = 1, De = 1.3252E3, m1 = 32.0419E-3, m2 = 39.948E-3, Re = 2.86E-10, vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3)
    y = teste.values()
    print(y[0], y[-1])
    plot(x[1:], y)
    show()
    '''
    teste = GibbsFreeEnergy(x, 1E-6, True)
    teste.setParameters(multiplicity = 1, De = 1.3252E3, m1 = 32.0419E-3, m2 = 39.948E-3, Re = 2.86E-10, vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3)
    y = teste.values()
    print(y[0], y[-1])
    plot(x[1:], y)
    show()
        
    
    
    