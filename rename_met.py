"""
Rename met_em files
"""

import glob
import os

files = glob.glob('met_em*')

for f in files:
    l = f.split(':')
    new = '%s_%s_%s' % tuple(l)
    os.system('mv %s %s' % (f, new))
