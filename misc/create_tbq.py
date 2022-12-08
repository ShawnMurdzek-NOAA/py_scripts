
import numpy as np

CQ=173.15-.05
R273=1./273.15
R61=6.1153*0.62198
ARP=77455.*41.9/461.525
BRP=64.*41.9/461.525

tbq = np.zeros(5002)
for k in range(1, 5002):
    CQ=CQ+.05
    evs=np.exp(17.67*(CQ-273.15)/(CQ-29.65))
    eis=np.exp(22.514-6.15E3/CQ)
    if CQ >= 273.15:
        tbq[k] = R61*evs
    else:
        tbq[k] = R61*eis
