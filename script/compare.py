#!/usr/bin/env python

import sys
import matplotlib
matplotlib.use('Agg')

from itertools import islice
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

DriverNet = {}
DriverNetPy = {}

# drivernet result
with open(sys.argv[1], 'r') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split()
        DriverNet[line[0]] = int(line[1])

# drivernetPy result
with open(sys.argv[2], 'r') as fin:
    for line in islice(fin, 1, None):
        line = line.strip().split('\t')
        DriverNetPy[line[1]] = int(line[0])


# venn
plt.figure(1)
plt.subplot(121)

venn2([set(list(DriverNet.keys())), set(list(DriverNetPy.keys()))], set_labels=('DriverNet', 'DriverNetPy'))
plt.title('Overlap')

# concordance
plt.subplot(122)
concordance_ratio = []

cumulative_count = 0
for ix, key in enumerate(DriverNetPy.keys()):
    if key in DriverNet:
        cumulative_count += 1
    concordance_ratio.append(cumulative_count/float(ix+1))

plt.plot(concordance_ratio)
plt.title('Concordance ratio')
plt.xlabel('Top ranked genes of DriverNetPy')
plt.ylabel('Concordance with DriverNet')

plt.tight_layout()
plt.savefig('compare.png', dpi = 300)
plt.clf()
plt.cla()
plt.close()

