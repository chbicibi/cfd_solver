import argparse
import glob
import os
import re
import shutil
import subprocess
from contextlib import contextmanager

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as plc

import myutils as ut
import mk_grid as mg
import cfd_main as cm

'''NOTE: u, v, p, f, w'''

i = 1
n = [2, 4][i]
nm = ['images_p', 'image_w'][i]
vr = [1, 0.1][i]
force = False
FFMPEG = next(filter(os.path.isfile, [
        r'C:\Tools\ffmpeg-4.2.2-win64-static\bin\ffmpeg.exe',
        os.path.abspath(r'bin\ffmpeg.exe')
    ]), None)


def plot(dest):
    flag = False
    with ut.chdir(dest):
        ut.mkdir(nm)
        # files = ut.iglobm('result/*.npy')
        files = ut.iglobm('*.dat')

        for file in files:
            basename = ut.basename(file, '.*')
            image = f'{nm}/{basename}.png'
            if os.path.isfile(image) and not force:
                continue
            print(dest, basename, end=' \r')

            # data = np.load(file)
            data = cm.read_raw(file)

            # print(data.shape)
            val = data[:, :, 3] # ((u, v, p, w), ...)
            if n == 2:
                a = val[val.shape[0]//2, 0]
                val -= a

            fig, ax = plt.subplots()
            fig.subplots_adjust(left=0.08, right=1, bottom=0, top=1)

            colors = [(0, '#ff0000'), (0.5, '#000000'), (1, '#00ff00')]
            cmap = plc.LinearSegmentedColormap.from_list('custom_cmap', colors)
            im = ax.imshow(val, cmap=cmap, vmin=-vr, vmax=vr)
            cax = fig.colorbar(im)

            ax_pos = ax.get_position()
            cax_pos0 = cax.ax.get_position()
            cax_pos1 = [cax_pos0.x0, ax_pos.y0, cax_pos0.x1 - cax_pos0.x0, ax_pos.y1 - ax_pos.y0]
            cax.ax.set_position(cax_pos1)

            # plt.show()
            fig.savefig(image)
            plt.close('all')
            flag = True
    return flag



def mk_v(dest, file='out.mp4'):
    if not FFMPEG:
        return
    with ut.chdir(os.path.join(dest, nm)):
        subprocess.call([FFMPEG, '-framerate', '30', '-y',
                        '-i', 'out_%04d.png', '-vcodec', 'libx264',
                        '-pix_fmt', 'yuv420p', '-r', '30', file])


def main():
    for d in os.listdir('.'):
        if os.path.isdir(d) and os.path.isfile(d+'/in2d.txt'):
            if plot(d):
                mk_v(d)


if __name__ == '__main__':
    main()

# set path=C:\Tools\ffmpeg-4.2.2-win64-static\bin;%path%
# ffmpeg -framerate 30 -i out_%04d.plt.png -vcodec libx264 -pix_fmt yuv420p -r 30 out.mp4
