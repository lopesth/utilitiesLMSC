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
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from pandas import read_csv
from random import randint
from numpy import exp
    
class CEPfitting(object):
    
    def __init__(self, function, csvFile, De, Re, xcol = 0, ycol = 1, xCorrection = 1, yCorrection = 1, method ='lm'):
        self.__method = method
        self.__De = De * yCorrection
        self.__Re = Re
        if function == "Improved Lennard-Jones":
            self.__guess = [randint(1, 2)]
            self.__func = self.LJ_Improved_ReqDe
        elif function == "Rydberg-6":
            self.__guess = [randint(1, 2) for i in range(6)]
            self.__func = self.rydberg6
        database = read_csv(csvFile)
        self.__yData = database.iloc[:,ycol] * yCorrection
        c = 0
        for y in self.__yData:
            if y < 0:
                break
            c += 1
        self.__yData = self.__yData[c-1:]
        self.__xData = database.iloc[:,xcol] * xCorrection
        self.__xData = self.__xData[c-1:]


            
    def LJ_Improved_ReqDe(self, x, beta):
        n_r = (beta + 4*((x/self.__Re)**2))
        return ((self.__De/(n_r - 6)) * ((6 * ((self.__Re/x)**(n_r))) - (n_r * ((self.__Re/x)**6))))

    def rydberg6(self, x, c1, c2, c3, c4, c5, c6):
        return -self.__De*(1 + c1*(x-self.__Re) + c2*(x-self.__Re)**2 + c3*(x-self.__Re)**3 + c4*(x-self.__Re)**4 + c5*(x-self.__Re)**5 + c6*(x-self.__Re)**6)*exp(-c1*(x-self.__Re))
        
    def setGuess(self, guess):
        self.__guess = guess
    
    def fit(self):
        self.__popt, self.__pcov = curve_fit(self.__func, self.__xData, self.__yData, p0=self.__guess, method=self.__method)
        
    def getOptimizedParameters(self):
        return self.__popt
    
    def getOptimizedValues(self):
        return self.__xData, self.__func(self.__xData, *self.__popt)

    def getReducedChi2(self):
        return chisquare(self.getOptimizedValues(), self.__yData)[0] / (len(self.__xData) - len(self.__popt))

    def getOriginalValues(self):
        return self.__xData, self.__yData

if __name__ == "__main__":
    
    fit1 = CEPfitting("Rydberg-6", '/Users/thiagolopes/Downloads/Ne_caso1.csv', De = 0.05652, Re = 2.60, xcol=0, ycol=7, yCorrection=43.3641)
    fit2 = CEPfitting("Rydberg-6", '/Users/thiagolopes/Downloads/He_caso2.csv', De = 0.05652, Re = 2.60, xcol=0, ycol=7, yCorrection=43.3641)
    fit1.fit() 
    fit2.fit()
    print(fit1.getOptimizedParameters())
    print(fit1.getReducedChi2())
    print(fit2.getOptimizedParameters())
    print(fit2.getReducedChi2())
    
