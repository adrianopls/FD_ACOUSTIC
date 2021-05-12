# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 15:46:45 2021

@author: Adriano Santana
"""

DIR_INPUT = "inputs"
DIR_OUTPUT = "outputs"

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import gridspec
from matplotlib.animation import FuncAnimation, writers

from calc_sg import staggeredGrid
from calc_rsg import rotatedStaggeredGrid
from animation import SGAnimation


np.seterr(all='raise')
#np.seterr(all='print')


#filename = 'B14_10_0600.rec.8bit.INVERTED-256px.png'
input_filename = 'B14_10_0600.rec.8bit-256px.png'
#
now = datetime.datetime.now()
output_filename = 'SG_' + now.strftime("%d%m%Y-%H%M%S") + '.mp4'
#

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
sim_steps = 1200
#
sou_x = 128
sou_y = 128
#
x_rec1 = 10
y_rec1 = 250
#
x_rec2 = 5
y_rec2 = 5
#
x_rec3 = 250
y_rec3 = 250
#
freq = 500000.0
amp = 1.0
#
dx = 15.0 #0.000296
dy = 15.0 #0.000296

#
dx = 0.000296
dy = 0.000296

#
# Read binary input file
input_vec = plt.imread(DIR_INPUT + '/' + input_filename)

# input_vec = np.zeros((256, 256))
# input_vec[100:156, :] = 1.0
# input_vec[65:80, :] = 1.0
# input_vec[40:50, :] = 1.0
# input_vec[176:190, :] = 1.0
# input_vec[205:215, :] = 1.0



# input_vec = np.zeros((256, 256))
# input_vec[0:50, :] = 1.0
# input_vec[100:156, :] = 1.0
# input_vec[200:256, :] = 1.0




# Calculate StaggeredGrid                                      0.000001
#wavefield, dx, dy, dt = rotatedStaggeredGrid(input_vec, sim_steps, dx=dx, dy=dy, f0=freq, amp=amp) #dx=0.000296, dy=0.000296)
wavefield, dx, dy, dt = staggeredGrid(input_vec, sim_steps, 
                                      sou_x=sou_x, sou_y=sou_y,
                                      dx=dx, dy=dy, f0=freq, amp=amp
                        ) #dx=0.000296, dy=0.000296)

#
sga_title = r'Staggered Grid Acoustic P-Wavefield $O(\Delta t^{2}, h^{4})$'
# Instantiating MPL animation base object
sga = SGAnimation(input_vec, wavefield, main_ax, ax0, ax1, ax2, dx=dx, dy=dy,
                  dt=dt, f0=freq, amp=amp, steps=sim_steps, 
                  sou_x=sou_x, sou_y=sou_y,
                  x_rec1=x_rec1, y_rec1=y_rec1,
                  x_rec2=x_rec2, y_rec2=y_rec2,
                  x_rec3=x_rec3, y_rec3=y_rec3,
                  title=sga_title)
# 
animation = FuncAnimation(fig, sga, frames=sim_steps, init_func=sga.init_func, 
                          interval=100, repeat=True, blit=False)

#
# print("\nSaving animation...")
# Saving output
#
# Based on example: https://www.youtube.com/watch?v=WXv7HQr_8SU
# 
Writer = writers['ffmpeg']
writer = Writer(fps=20, metadata={'artist': 'Adriano Santana', 'year': '2021'},
                bitrate=1800)

animation.save(DIR_OUTPUT + '/' + output_filename, writer)
print("Saved to " + DIR_OUTPUT + '/' + output_filename)

# plt.show()



