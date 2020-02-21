#! /usr/bin/env python3

import argparse
import glob
import os
import re
import shutil
import subprocess
from contextlib import contextmanager
import numpy as np
import myutils as ut
import mk_grid as mg

'''NOTE: u, v, p, f, w'''

EXE_EXT = 'exe' if os.environ.get('OS') == 'Windows_NT' else 'out'

def run_cfd(resume=False):
    info = mg.get_info('in2d.txt')
    dest = info['dest']
    os.makedirs(dest, exist_ok=True)

    if resume:
        if not rewrite_resume(dest, f'in2d.txt'):
            return

    for file in ('in2d.txt', 'grid.csv', 'grid.png'):
        shutil.copy(file, dest)

    exe = glob.glob(f'bin/a.{EXE_EXT}')[0]
    res = subprocess.run(['./'+ exe])
    print(res.returncode)
    return res.returncode


################################################################################

def rewrite_resume(dest, file):
    with ut.chdir(dest):
        plts = glob.glob('*.plt')

    if not plts:
        return False

    last = max(map(int, (re.search(r'\d+', f)[0] for f in plts)))

    with open(file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    iout = int(re.findall(r'\d+', lines[30])[0])

    with open(file, 'w', encoding='utf-8') as f:
        for i, l in enumerate(lines):
            if i == 15:
                ll = re.split(r'(?<= )(?=\S)', l)
                ll[1] = f'{(last+1)*iout:<9} '
                l = ''.join(ll)

            elif i == 27:
                ll = re.split(r'(?<= )(?=\S)', l)
                ll[4] = f'{last+1}\n'
                l = ''.join(ll)

            f.write(l)
    return True


def read_plt(file):
    prog = re.compile(r'\d+')
    print('reading:', file, end='\r')
    with open(file, 'r') as f:
        f.readline()
        line = f.readline()
        nx, ny = map(int, prog.findall(line))
        res = np.empty((ny, nx, 5), np.float32)

        for j in range(ny):
          for i in range(nx):
            line = f.readline()
            res[j, i, :] = [float(s) for s in line.split()][2:]
    return res


@contextmanager
def post_base(dest=None):
    if not dest:
        info = mg.get_info('in2d.txt')
        dest = info['dest']
    with ut.chdir(dest):
        yield


def collect_result(dest):
    print('collect_result')
    odir = 'result'
    rdir = '__raw__'
    with post_base(dest):
        os.makedirs(odir, exist_ok=True)
        os.makedirs(rdir, exist_ok=True)
        files = sorted(glob.iglob('out_*.plt'))
        for i, file in enumerate(files):
            # print(i, end='\r')
            ofile = os.path.join(odir, os.path.basename(file) + '.npy')
            if os.path.isfile(ofile):
                continue
            data = read_plt(file)
            np.save(ofile, data)
            shutil.move(file, rdir)


def pack_data():
    def f_(a):
        return a[:, :, 4] # => vorticity
    odir = 'result'
    with post_base():
        with ut.chdir(odir):
            files = sorted(glob.iglob('out_*.npy'))
            # data = [np.load(file) for file in files]
            data = [f_(np.load(file)) for file in files[0:1000:10]]
        # ofile = f'vorticity_{len(files)}.npy'
        ofile = 'vorticity_100.npy'
        print(ofile)
        np.save(ofile, data)

################################################################################

def __test__():
    '''+-0.5~1.0'''
    from PIL import Image

    def f_(a_mono):
        r = np.clip(-5 * 256 * a_mono, 0, 255)
        g = np.clip(5 * 256 * a_mono, 0, 255)
        b = np.zeros_like(a_mono, dtype=np.uint8)
        a = [c.astype(dtype=np.uint8)[:, :, np.newaxis] for c in (r, g, b)]
        img = np.concatenate(a, axis=2)
        return img

    # file = glob.glob('plate_20/all_data*.npy')[0]
    # for file in glob.glob('out_*.npy'):
    for file in ('vorticity_100_a0.npy',):
        ofile = file + '.png'
        # if os.path.exists(ofile):
        #     continue
        print(file)
        data = np.load(file)
        img = f_(data[-1])
        im = Image.fromarray(img)
        # print(data.shape)
        # im.show()
        im.save(ofile)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-run', action='store_true', help='start calculating')
    parser.add_argument('-res', action='store_true', help='converting result files')
    parser.add_argument('--resume', action='store_true',
                        help='resume')
    parser.add_argument('-dest', help='output directory name')
    parser.add_argument('-pack', action='store_true', help='packing result files')
    parser.add_argument('-test', action='store_true', help='test mode')
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    if args.run:
        with ut.stopwatch('CFD'):
            run_cfd(resume=args.resume)
    elif args.res:
        with ut.stopwatch('RES'):
            if args.dest:
                collect_result(args.dest)
            else:
                for d in os.listdir('.'):
                    if os.path.isdir(d) and os.path.isfile(d+'/in2d.txt'):
                        print(d)
                        collect_result(d)
    elif args.pack:
        pack_data()
    elif args.test:
        __test__()


if __name__ == '__main__':
    main()
