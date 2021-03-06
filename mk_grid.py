import argparse
import csv
import math
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from PIL import Image, ImageFilter
from scipy import interpolate

from plot_grid import plot_shape


################################################################################

def rotate(x, y, cx, cy, a):
    x_ = (x - cx) * math.cos(a) - (y - cy) * math.sin(a) + cx
    y_ = (x - cx) * math.sin(a) + (y - cy) * math.cos(a) + cy
    return x_, y_


def naca4(mr, pr, t, c, n=1000, get_func=False):
    m = mr * c # 最大キャンバー
    p = pr * c # 最大キャンバー位置
    ct5 = 5 * c * t

    def yc(x):
        if x < p:
            return m * x       / p       ** 2 * ( 2 * p -     x / c)
        else:
            return m * (c - x) / (1 - p) ** 2 * (-2 * p + 1 + x / c)

    def dyc_dx(x):
        if x < p:
            return 2 * m / p       ** 2 * (p - x / c)
        else:
            return 2 * m / (1 - p) ** 2 * (p - x / c)

    def yt(x):
        xc = x / c
        return ct5 * (  0.2969 * xc ** 0.5
                      - 0.1260 * xc
                      - 0.3516 * xc ** 2
                      + 0.2843 * xc ** 3
                      - 0.1015 * xc ** 4)

    def shape(x):
        ytx = yt(x)
        ycx = yc(x)
        tt = math.atan(dyc_dx(x))
        yts = ytx * math.sin(tt)
        ytc = ytx * math.cos(tt)
        xu, xl, yu, yl = x - yts, x + yts, ycx + ytc, ycx - ytc
        return xu, xl, yu, yl

    def line_f(xus, xls, yus, yls):
        line_u = interpolate.interp1d(xus, yus, kind='cubic')
        line_l = interpolate.interp1d(xls, yls, kind='cubic')
        xlim = max(xus[0], xls[0]), min(xus[-1], xls[-1])

        # xs = np.linspace(*xlim, 1000)
        # print(xs)
        # yu = line_u(xs)
        # yl = line_l(xs)
        # print(xus[-1], yus[-1])
        # print(xls[-1], yls[-1])
        # plt.plot(xs, yu)
        # plt.plot(xs, yl)
        # plt.show()

        def f_(x, y):
            # print(x, y)
            if xlim[0] <= x <= xlim[1]:
                return line_l(x) <= y <= line_u(x)
            else:
                return False
        return f_

    if get_func:
        return line_f(*zip(*map(shape, np.linspace(0, 1, n))))
    else:
        return list(map(shape, np.linspace(0, 1, n)))


################################################################################

def grid_plain(shape):
    grid = np.ones(shape, dtype=np.uint8)
    return grid


def grid_plate(shape):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)
    grid[ny//2-2:ny//2+2, nx//8:nx//8+ny] = 0
    return grid


def grid_plate_a(shape, alpha):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    left = nx // 8
    length = nx // 4
    thickness = 4

    a_rad = math.radians(alpha)
    cx = left + length // 2
    cy = ny // 2
    base_x = left, length
    base_y = cy - thickness // 2, thickness

    for i in range(1000): # x
        for j in range(1000): # y
            x = base_x[0] + base_x[1] * i / 999
            y = base_y[0] + base_y[1] * j / 999
            r = rotate(x, y, cx, cy, a_rad)
            x_, y_ = map(int, r)
            grid[y_, x_] = 0
    return grid


def grid_wing(shape, alpha, params=(0.02, 0.4, 0.12, 1.0)):
    if type(params) is str:
        params = (int(params[0]) * 0.01,
                  int(params[1]) * 0.1,
                  int(params[2:]) * 0.01,
                  1.0)

    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    left = nx // 8
    length = nx // 4

    a_rad = math.radians(alpha)
    cx = left + length // 2
    cy = ny // 2

    # for xu, xl, yu, yl in naca4(*params, n=1000):
    #     for x0, y0 in ((xu, yu), (xl, yl)):
    #         for i in np.linspace(0, 1, 100):
    #             x = left + length * x0
    #             y = cy - length * y0 * i
    #             x_, y_ = map(int, rotate(x, y, cx, cy, a_rad))
    #             grid[y_, x_] = 0

    naca_f = naca4(*params, n=1000, get_func=True)
    for i in range(nx):
        for j in range(ny):
            x_, y_ = rotate(i+0.5, ny-j+0.5, cx, cy, a_rad)
            ix = (x_ - left) / length
            iy = (y_ - cy) / length

            if naca_f(ix, iy):
                grid[j, i] = 0

    return grid


################################################################################

def grid_cylinder(shape, alpha=0, ry=6):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    r = ny / ry
    cx = nx / 4
    cy = ny / 2

    for i in range(nx):
        for j in range(ny):
            if (cx - i - 0.5) ** 2 + (cy - j - 0.5) ** 2 < r ** 2:
                grid[j, i] = 0
    return grid


def grid_ellipse(shape, alpha=0, ry=8):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    a = ny / 6
    b = ny / ry
    cx = nx / 4
    cy = ny / 2
    a_rad = math.radians(-alpha)

    for i in range(nx):
        for j in range(ny):
            x_, y_ = rotate(i + 0.5, j + 0.5, cx, cy, a_rad)
            if ((cx - x_) / a) ** 2 + ((cy - y_) / b) ** 2 < 1:
                grid[j, i] = 0
    return grid


def grid_prism(shape, alpha=0, ry=8):
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    w = ny / 6 / math.sqrt(2) # half
    h = w * 6 / ry
    cx = nx / 4
    cy = ny / 2
    a_rad = math.radians(-alpha)

    for i in range(nx):
        for j in range(ny):
            x_, y_ = rotate(i + 0.5, j + 0.5, cx, cy, a_rad)
            if cx - w < x_ < cx + w and cy - h < y_ < cy + h:
                grid[j, i] = 0
    return grid


################################################################################

def get_info(input_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    nx, ny = [int(s) for s in re.findall(r'\d+', lines[24])][:2]
    dest = re.findall(r'\S+', lines[30])[1].strip()
    info = {
        'nx': nx,
        'ny': ny,
        'dest': dest
    }
    return info


def dump_grid(grid, file):
    with open(file, 'w', encoding='utf-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerows(grid)


################################################################################

def make_grid(args, info):
    a = args.a
    t = args.t
    print('alpha=', a, 'deg', 'type=', t)
    shape = [info[k] for k in ('ny', 'nx')]

    if t == 'wing':
        grid = grid_wing(shape, a, args.naca)

    elif t == 'cylinder':
        grid = grid_cylinder(shape, a, ry=args.r)

    elif t == 'ellipse':
        grid = grid_ellipse(shape, a, ry=args.r)

    elif t == 'prism':
        grid = grid_prism(shape, a, ry=args.r)

    elif t == 'plate':
        grid = grid_plate_a(shape, a)

    else:
        raise AttributeError

    return grid


################################################################################

def cvt_grid(file):
    if not file.endswith('.png'):
        file += '.png'

    if not os.path.isfile(file):
        print('pngファイルがありません:', file)
        raise FileNotFoundError

    image = Image.open(file)
    array = np.asarray(image)

    if array.ndim > 2:
        array = array[:, :, 0]

    print(array.shape, array.dtype)

    grid = np.where(array < 128, 1, 0).astype('u1')

    return grid


################################################################################

def get_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument('live_id', help='live id', nargs='*')
    # parser.add_argument('-f', help='get live ids from file', action='store_true')
    # parser.add_argument('-w', help='watch live page', action='store_true')
    # parser.add_argument('-c', help='getting comments', action='store_true')
    # parser.add_argument('-dl', help='call main', action='store_true')
    # parser.add_argument('-ck', help='check filename', action='store_true')
    # parser.add_argument('-op', help='open live page', action='store_true')
    # parser.add_argument('-debug', help='enable debug mode', action='store_true')
    # parser.add_argument('-concat', help='concatenate video files', action='store_true')

    parser.add_argument('-a', default=0, type=float, help='alpha')
    parser.add_argument('-t', choices=['wing', 'cylinder', 'plate', 'ellipse',
                                       'prism'], default='cylinder',
                        help='grid shape')
    parser.add_argument('-r', default=6, type=float, help='y ratio')
    parser.add_argument('-i', default='', help='create from png')
    parser.add_argument('--naca', default='0012', help='params for NACA4')
    parser.add_argument('--plot', action='store_true', help='plot shape')
    parser.add_argument('--test', action='store_true', help='test mode')
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    if args.plot:
        plot_shape('grid.csv', 'grid.png')
        return

    info = get_info('in2d.txt')

    if args.i:
        grid = cvt_grid(args.i)
    else:
        grid = make_grid(args, info)

    print(grid.shape, grid.dtype)
    dump_grid(grid, 'grid.csv')
    plot_shape('grid.csv', 'grid.png')


if __name__ == '__main__':
    main()
