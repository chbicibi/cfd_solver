#! /usr/bin/env python3

import argparse
import csv
import math
import re
import numpy as np


################################################################################

def rotate(x, y, cx, cy, a):
    x_ = (x - cx) * math.cos(a) - (y - cy) * math.sin(a) + cx
    y_ = (x - cx) * math.sin(a) + (y - cy) * math.cos(a) + cy
    return x_, y_


def naca4(mr, pr, t, c, n=1000):
    m = mr * c
    p = pr * c
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
        return x, ycx, x - yts, x + yts, ycx + ytc, ycx - ytc
    return [shape(i / (n - 1.0)) for i in range(n)]


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
    ny, nx = shape
    grid = np.ones(shape, dtype=np.uint8)

    left = nx // 8
    length = nx // 4

    a_rad = math.radians(alpha)
    cx = left + length // 2
    cy = ny // 2

    for _x, _yc, xu, xl, yu, yl in naca4(*params, n=1000):
        for x0, y0 in ((xu, yu), (xl, yl)):
            for i in range(100):
                x = left + length * x0
                y = cy - length * y0 * i / 99
                x_, y_ = map(int, rotate(x, y, cx, cy, a_rad))
                grid[y_, x_] = 0
    return grid


################################################################################

def get_info(input_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    nx, ny = [int(s) for s in re.findall(r'\d+', lines[24])]
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
        grid = grid_wing(shape, a)
    else:
        grid = grid_plate_a(shape, a)
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

    parser.add_argument('-a', default=0, type=int, help='alpha')
    parser.add_argument('-t', default='plate', type=str, help='grid type')
    parser.add_argument('-test', action='store_true', help='test mode')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    info = get_info('in2d.txt')
    grid = make_grid(args, info)
    print(grid.shape)
    dump_grid(grid, 'grid.csv')


if __name__ == '__main__':
    main()
