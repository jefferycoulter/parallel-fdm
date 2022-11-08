#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# dt = 0.05

# ======================== create video of simulation ========================
# print("loading data")
# data = pd.read_csv("../data/data.csv", dtype=float, delimiter=",")
# print("converting to numpy array")
# data = data.to_numpy()
# print("deleting column")
# data = np.delete(data, -1, 1)
# print("reshaping array")
# data = data.reshape(5000, 100, 100)

# fig_a, ax_a = plt.subplots()
# cbar_ax = fig_a.add_axes([0.80, 0.15, 0.05, 0.70])
# ims = []
# print("creating plots")
# for i in range(len(data)):
#     im = ax_a.imshow(data[i][:][:], origin="lower", cmap="hot", animated=True)
#     fig_a.colorbar(im, cbar_ax)
#     if i == 0:
#         im = ax_a.imshow(data[i][:][:], origin="lower", cmap="hot", animated=True)
#     ims.append([im])

# ani = animation.ArtistAnimation(fig_a, ims, interval=5, blit=True, repeat_delay=1000)
# ani.save("../graphics/animations/movie.gif")
        

# ============================ plot shape array 2D ===========================

# shape = np.genfromtxt("../data/shape.csv", dtype=float, delimiter=",")
# #shape = shape.to_numpy()
# shape = np.delete(shape, -1, 0)
# shape = shape.reshape(100, 100)

# fig, ax = plt.subplots()
# fig.tight_layout()
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
# im = ax.imshow(shape, origin="lower", cmap="binary")
# fig.colorbar(im, cbar_ax)

# ============================ plot shape array 3D ===========================
print("loading u data")
u = np.genfromtxt("../data/data_u.csv", dtype=float, delimiter=",")
u = np.delete(u,  -1, 0)
u = np.delete(u,  -1, 1)
u = u.reshape(40, 100, 100, 100)

print("loading v data")
v = np.genfromtxt("../data/data_v.csv", dtype=float, delimiter=",")
v = np.delete(v,  -1, 0)
v = np.delete(v,  -1, 1)
v = v.reshape(40, 100, 100, 100)

print("loading uv data")
uv = np.genfromtxt("../data/data_uv.csv", dtype=float, delimiter=",")
uv = np.delete(uv,  -1, 0)
uv = np.delete(uv,  -1, 1)
uv = uv.reshape(40, 100, 100, 100)

fig, ax = plt.subplots()
fig.tight_layout()
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
im0 = ax.imshow(u[10][50][:], origin="lower", cmap="Blues")
im1 = ax.imshow(v[10][50][:], alpha=0.5, origin="lower", cmap="Greens")
im2 = ax.imshow(uv[10][50][:], alpha=0.5, origin="lower", cmap="Purples")
#fig.colorbar(im, cbar_ax)
plt.show()

# shape = np.genfromtxt("../data/shape.csv", dtype=float, delimiter=",")
# shape = np.delete(shape, -1, 0)
# shape = shape.reshape(100, 100, 100)

# fig, ax = plt.subplots()
# fig.tight_layout()
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.70])
# im = ax.imshow(shape[50][:], origin="lower", cmap="binary")
# fig.colorbar(im, cbar_ax)

# shape = np.genfromtxt("../data/shape.csv", dtype=float, delimiter=",")
# #shape = shape.to_numpy()
# shape = np.delete(shape, -1, 0)
# shape = shape.reshape(100, 100, 100)

# fig = plt.subplots()
# ax = plt.axes(projection="3d")
# ax.voxels(shape)

# ============================== create video 3D =============================

# print("loading data")
# shape = np.genfromtxt("../data/data.csv", dtype=float, delimiter=",")
# print("deleting column")
# data = np.delete(data, -1, 1)
# print("reshaping array")
# data = data.reshape(1000, 100, 100, 100)

# fig_a = plt.subplots()
# ax = plt.axes(projection="3d")
# ims = []
# print("creating plots")
# for i in range(len(data)):
#     im = ax.vocels(data[i][:][:][:], origin="lower", cmap="hot", animated=True)
#     if i == 0:
#         im = ax.vocels(data[i][:][:][:], origin="lower", cmap="hot", animated=True)
#     ims.append([im])

# ani = animation.ArtistAnimation(fig_a, ims, interval=5, blit=True, repeat_delay=1000)
# ani.save("../graphics/animations/movie.gif")


