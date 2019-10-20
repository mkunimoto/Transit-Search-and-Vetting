import sys
import numpy as np

g1 = float(sys.argv[1])
g2 = float(sys.argv[2])
n = len(sys.argv) - 2
p = n/3
per = sys.argv[3:3+p]
per = np.array([float(i) for i in per])
epo = sys.argv[3+p:3+p*2]
epo = np.array([float(i) for i in epo])
dep = sys.argv[3+p*2:3+3*p]
dep = np.array([float(i) for i in dep])
rat = np.sqrt(dep)

f = open('LMfit_template.dat','w')

f.write('RHO 1.408 0.0 0.0 0.0' + '\n') # assumes solar density
f.write('NL1 ' + str(g1) + ' 0.0 0.0 0.0' + '\n')
f.write('NL2 ' + str(g2) + ' 0.0 0.0 0.0' + '\n')
f.write('NL3 0.0 0.0 0.0 0.0' + '\n')
f.write('NL4 0.0 0.0 0.0 0.0' + '\n')
f.write('DIL 0.0 0.0 0.0 0.0' + '\n')
f.write('VOF 0.0 0.0 0.0 0.0' + '\n')
f.write('ZPT 0.0 0.0 0.0 0.0' + '\n')
for i in range(p):
    f.write('EP1 ' + str(epo[i]) + ' 0.0 0.0 0.0' + '\n')
    f.write('PE1 ' + str(per[i]) + ' 0.0 0.0 0.0' + '\n')
    f.write('BB1 0.0 0.0 0.0 0.0' + '\n')
    f.write('RD1 ' + str(rat[i]) + ' 0.0 0.0 0.0' + '\n')
    f.write('EC1 0.0 0.0 0.0 0.0' + '\n')
    f.write('ES1 0.0 0.0 0.0 0.0' + '\n')
    f.write('KR1 0.0 0.0 0.0 0.0' + '\n')
    f.write('TE1 0.0 0.0 0.0 0.0' + '\n')
    f.write('EL1 0.0 0.0 0.0 0.0' + '\n')
    f.write('AL1 0.0 0.0 0.0 0.0' + '\n')

f.close()
