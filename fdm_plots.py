#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt

dt = 0.05
data = np.genfromtxt("data/data.csv", dtype=float, delimiter=",")

data = np.delete(data, -1, 1)
data = data.reshape(5001, 200, 200)
        

fig, ax = plt.subplots(nrows=2, ncols=2)
fig.tight_layout()

ax[0,0].set_title("time = {:.2f}".format(dt*10))
ax[0,0].imshow(data[10][:][:])
ax[0,1].set_title("time = {:.2f}".format(dt*1000))
ax[0,1].imshow(data[250][:][:])
ax[1,0].set_title("time = {:.2f}".format(dt*2500))
ax[1,0].imshow(data[500][:][:])
ax[1,1].set_title("time = {:.2f}".format(dt*5000))
ax[1,1].imshow(data[1000][:][:])

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(ax[1,1], cax=cbar_ax)
