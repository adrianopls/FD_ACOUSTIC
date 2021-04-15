# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 15:46:45 2021

@author: Adriano Santana
"""

DIR_INPUT = "inputs"
DIR_OUTPUT = "outputs"


import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.animation import FuncAnimation, writers

from calc import staggeredGrid
from animation import SGAnimation


#filename = 'B14_10_0600.rec.8bit.INVERTED-256px.png'
input_filename = 'B14_10_0600.rec.8bit-256px.png'
#
output_filename = 'simula.mp4'


# Initialize the plot
fig = plt.figure(figsize=(16,12))
gs = gridspec.GridSpec(ncols=6, nrows=6, figure=fig)
#
main_ax = fig.add_subplot(gs[:, 0:4])
#
ax0 = fig.add_subplot(gs[0:2, 4:6])
ax1 = fig.add_subplot(gs[2:4, 4:6])
ax2 = fig.add_subplot(gs[4:6, 4:6])
#
sim_steps = 60
sou_x = 128
sou_y = 128
# Read binary input file
input_vec = plt.imread(DIR_INPUT + '/' + input_filename)
# Calculate StaggeredGrid 
wavefield, dx, dy, dt = staggeredGrid(input_vec, sim_steps)
# Instantiating MPL animation base object
sga = SGAnimation(input_vec, wavefield, main_ax, ax0, ax1, ax2, dx=dx, dy=dy,
                  dt=dt, f0=500000, amp=30.0, steps=sim_steps, 
                  sou_x=sou_x, sou_y=sou_y)
# 
animation = FuncAnimation(fig, sga, frames=sim_steps, init_func=sga.init_func, 
                          interval=50, repeat=False, blit=False)
# Saving output
#
# Based on example: https://www.youtube.com/watch?v=WXv7HQr_8SU
# 
# Writer = writers['ffmpeg']
# writer = Writer(fps=20, metadata={'artist': 'Adriano Santana', 'year': '2021'},
#                 bitrate=1800)
# #
# animation.save(DIR_OUTPUT + '/' + output_filename, writer)
#
plt.show()




PAREI NA WAVELET.SHAPE ESTRANHA: (33783784,)


