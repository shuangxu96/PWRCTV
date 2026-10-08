"""
Python based PWRCTV.
"""
import os
import sys
import time

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src import pwrctv, warm_up
from utils.make_noisy import make_noisy

print("warm-up (JIT compile) ...")
t0 = time.perf_counter()
warm_up(verbose=True)
dt = time.perf_counter() - t0
print(f'warm up %f seconds'%(dt))

for dn in ("Florence", "Milan"):
    
    for case in (1, 2, 3, 4, 5):
        Nhsi, Ohsi, Pan = make_noisy(dn, case)

        tt = 0.4 if case == 1 else 0.7
        qq = 10 if case == 1 else 5

        t0 = time.perf_counter()
        out, _, _ = pwrctv(Nhsi, Pan, 100, 1, np.array([tt, tt]), 4, qq)
        dt = time.perf_counter() - t0
        
        print(dn, ' case',case, dt, 'seconds')

