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

def run_cfd(exe_path, resume=''):
    if not os.path.isfile(exe_path):
        raise FileNotFoundError(exe_path)

    if resume:
        if not rewrite_resume(resume):
            return
    else:
        info = mg.get_info('in2d.txt')
        dest = info['dest']
        os.makedirs(dest, exist_ok=True)

        for file in ('in2d.txt', 'grid.csv', 'grid.png'):
            shutil.copy(file, dest)

    # exe = glob.glob(exe_path)[0]
    res = subprocess.run(['./'+ exe_path])
    print(res.returncode)
    return res.returncode


################################################################################

def rewrite_resume(resume):
    if not os.path.isdir(resume):
        raise FileNotFoundError(resume)

    with ut.chdir(resume):
        plts = glob.glob('*.plt')

        if not plts:
            return False

        last = max(map(int, (re.search(r'\d+', f)[0] for f in plts)))

        with open('in2d.txt', 'r', encoding='utf-8') as f:
            lines = f.readlines()
        iout = int(re.findall(r'\d+', lines[30])[0])

        with open('in2d.txt', 'w', encoding='utf-8') as f:
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

        for file in ('in2d.txt', 'grid.csv', 'grid.png'):
            shutil.copy(file, '../')

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


def read_raw(file):
    def readf():
        with open(file, 'rb') as f:
            nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
            u = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))
            yield (u[1:-1, :-1] + u[1:-1, 1:]) / 2

            nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
            v = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))
            yield (v[:-1, 1:-1] + v[1:, 1:-1]) / 2

            nx, ny, nb = (int.from_bytes(f.read(4), 'little') for i in range(3))
            p = np.frombuffer(f.read(nx*ny*nb), dtype=f'f{nb}').reshape((ny, nx))
            yield p[1:-1, 1:-1]

            w = (u[:-2, :-1] + u[:-2, 1:] - u[2:, :-1] - u[2:, 1:] +
                 v[:-1, 2:] + v[1:, 2:] - v[1:, :-2] - v[:-1, :-2]) / 4
            yield w

    return np.stack(list(readf())).transpose(1, 2, 0) # (C, H, W) -> (H, W, C)


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
        files = ut.fsort(glob.iglob('out_*.plt'))
        last = len(files) - 1
        if last < 1:
            return
        os.makedirs(odir, exist_ok=True)
        os.makedirs(rdir, exist_ok=True)
        for i, file in enumerate(files):
            # print(i, end='\r')
            ofile = os.path.join(odir, os.path.basename(file) + '.npy')
            if not os.path.isfile(ofile):
                data = read_plt(file)
                np.save(ofile, data)
            if i < last:
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
    parser.add_argument('-exe', default='a', help='start calculating')
    parser.add_argument('-res', action='store_true', help='converting result files')
    parser.add_argument('--resume', '-r', default='',
                        help='resume directory name')
    parser.add_argument('-dest', help='output directory name')
    parser.add_argument('-pack', action='store_true', help='packing result files')
    parser.add_argument('-test', action='store_true', help='test mode')
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    if args.run:
        exe_path = f'bin/{args.exe}'
        exe_ext = '.exe' if os.environ.get('OS') == 'Windows_NT' else '.out'
        if not exe_path.endswith(exe_ext):
            exe_path += exe_ext
        with ut.stopwatch('CFD'):
            run_cfd(exe_path, resume=args.resume)
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
