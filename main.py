from APP.termo import *
from APP.cepFit import *
import numpy as np
import matplotlib.pyplot as pyplot
import re, os


padraoCSV = re.compile('[\W\w]+.csv$')
arqv = [x for x in os.listdir() if padraoCSV.match(x)!=None]
print(arqv)
for arq in arqv:
    ycol = 7
    label = 'SAPT2/aug-cc-pVTZ'
    if 'He' in arq:
        if 'caso1' in arq:
            Re = 3.25
            De = 0.08952063
        elif 'caso2' in arq:
            Re = 2.60
            De = 0.05651687
        else:
            Re = 3.20
            De = 0.09353856      
    elif 'Ne' in arq:
        if 'caso1' in arq:
            Re = 3.35
            De = 0.15999027
        elif 'caso2' in arq:
            Re = 2.6
            De = 0.12909926
        else:
            Re = 3.25
            De = 0.18134327
    elif 'Ar' in arq:
        label = 'SAPT2/aug-cc-pVDZ'
        if 'caso1' in arq:
            Re = 3.70
            De = 0.36666105
        elif 'caso2' in arq:
            Re = 2.85
            De = 0.31672613
        else:
            Re = 3.55
            De = 0.42878623 
    else:
        label = 'SAPT2+(3)/aug-cc-pVDZ'
        ycol = 9
        if 'caso1' in arq:
            Re = 3.85
            De = 0.42768406
        elif 'caso2' in arq:
            Re = 2.95
            De = 0.40148747
        else:
            Re = 3.70
            De = 0.48679085
    fit_LJ = CEPfitting("Improved Lennard-Jones", arq, De, Re, ycol=ycol, yCorrection=43.3641)
    fit_LJ.fit()
    fit_R6 = CEPfitting("Rydberg-6", arq, De, Re, ycol=ycol, yCorrection=43.3641)
    fit_R6.fit()

    pyplot.subplot(111, label=arq)
    pyplot.plot(*fit_LJ.getOriginalValues(), 'o', color='#000000', markersize=4, label=label)
    pyplot.plot(*fit_LJ.getOptimizedValues(), '-', color='#0000cc', markersize=4, label='ILJ')
    pyplot.plot(*fit_R6.getOptimizedValues(), '-', color='#cc0000', markersize=4, label='R6')
    pyplot.legend(loc='best', frameon=False)
    pyplot.tight_layout()
    pyplot.savefig(arq.replace('.csv', '.png'), dpi=500, transparent=True)         

'''
Temp = np.arange(90, 200, 1)

Sis = GibbsFreeEnergy(Temp, 1E-6, True)
Sis.setParameters(multiplicity = 1, De = 1.3252E3, m1 = 32.0419E-3, m2 = 39.948E-3, Re = 2.86E-10, vMax = 7, w_e = 0.40307003819632E3, w_exe = 0.035017694E3, w_eye = 0.0004824939203489E3)
y = Sis.values()
y = np.array(y)
x = Temp[1:]


pyplot.figure(1)
pyplot.plot(x, y*0.001, '-', color='#0000cc')
pyplot.show()'''