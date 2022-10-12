#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dt = 0.05
# data = np.genfromtxt("../data/data.csv", dtype=float, delimiter=",")
print("loading data")
data = pd.read_csv("../data/data.csv", dtype=float, delimiter=",")
print("converting to numpy array")
data = data.to_numpy()
print("deleting column")
data = np.delete(data, -1, 1)
print("reshaping array")
data = data.reshape(5000, 200, 200)
        
#fig, ax = plt.subplots()
# =============================================================================

# shape = np.genfromtxt("../data/shape.csv", dtype=float, delimiter=",")
# shape = np.delete(shape, -1, 0)
# shape = shape.reshape(200,200)
# ax.imshow(shape, origin="lower", cmap="binary")
# =============================================================================



# ax[0,0].set_title("time = {:.2f}".format(dt*0))
# im = ax[0,0].imshow(data[0][:][:], cmap="hot")
# ax[0,1].set_title("time = {:.2f}".format(dt*1000))
# ax[0,1].imshow(data[250][:][:], cmap="hot")
# ax[1,0].set_title("time = {:.2f}".format(dt*2500))
# ax[1,0].imshow(data[500][:][:], cmap="hot")
# ax[1,1].set_title("time = {:.2f}".format(dt*5000))
# ax[1,1].imshow(data[1000][:][:], cmap="hot")

# fig.tight_layout()
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# fig.colorbar(im, cax=cbar_ax)


fig_a, ax_a = plt.subplots()
cbar_ax = fig_a.add_axes([0.80, 0.15, 0.05, 0.70])
ims = []
print("creating plots")
for i in range(len(data)):
    im = ax_a.imshow(data[i][:][:], origin="lower", cmap="hot", animated=True)
    fig_a.colorbar(im, cbar_ax)
    if i == 0:
        im = ax_a.imshow(data[i][:][:], origin="lower", cmap="hot", animated=True)
    ims.append([im])

ani = animation.ArtistAnimation(fig_a, ims, interval=5, blit=True, repeat_delay=1000)
ani.save("../animations/movie.gif")

