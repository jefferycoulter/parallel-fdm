#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def ComputeSpeedup(f, p):
    return 1 / ((1 - f) + (f / p))

# ratio of code that is parallelizable to the entire program runtime
f = 653801 / 700009 # [us]

times = [9, 5, 4, 4, 9, 24]

nproc = [1, 2, 4, 6, 8, 10]

fig, ax = plt.subplots()
fig.suptitle("Execution time vs number of processors")
ax.set_title("(grid size 100x100x100)")
ax.set_xlabel("number of processors")
ax.set_ylabel("time [s]")
ax.plot(nproc, times)


