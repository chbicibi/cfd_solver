import argparse
import math
import os
import re
import subprocess
import time
from tqdm import tqdm
from zipfile import ZipFile, ZIP_DEFLATED
import numpy as np
import myutils as ut


def comp_main():
    def read(files):
        for file in files:
            name = re.search(r'out_\d+', file)[0]
            print(name, end=' \r')
            yield name, np.load(file)

    for d in ut.fsort(os.listdir('.')):
        loc = f'{d}/result'
        if not os.path.isdir(loc):
            continue

        try:
            with ut.chdir(loc):
                # if any(ut.iglobm('*.npz')):
                #     continue

                if any(ut.iglobm('*.npy')):
                    print(loc)

                for i in range(10):
                    out = f'out_{i:03d}.npz'
                    if os.path.isfile(out):
                        continue

                    l = i * 100
                    u = (i + 1) * 100 if i < 9 else 1001
                    files = [f'out_{j:04d}.plt.npy' for j in range(l, u)]

                    if not all(map(os.path.isfile, files)):
                        continue

                    print(f'{l}=>{u}', ' '*20, end='\r')
                    h = dict(read(files))

                    try:
                        np.savez_compressed(out, **h)

                    except:
                        if os.path.isfile(out):
                            os.remove(out)
                        raise

                    ck = np.load(out)
                    for nm in ck.files:
                        f = f'{nm}.plt.npy'
                        if os.path.isfile(f):
                            os.remove(f)

        except OSError:
            print('OSError')
            time.sleep(10)
            continue


def zip_main():
    for loc in ut.fsort(ut.iglobm('**/image_*')):
        # loc = f'{d}/image_w'
        if not os.path.isdir(loc):
            continue

        with ut.chdir(loc):
            file = f'{ut.basename(loc)}.zip'
            if os.path.isfile(file):
                continue

            print(loc)
            # subprocess.run('zip image_w *.png', shell=True)

            try:
                def f_():
                    for png in ut.fsort(ut.iglobm('*.png')):
                        print(png, os.path.getsize(png), end='\r')
                        with open(png, 'rb'):
                            yield png

                pngs = list(f_())
                if not pngs:
                    continue

                with ZipFile(file, 'w', compression=ZIP_DEFLATED) as z:
                    for png in pngs:
                        z.write(png)
            except:
                if os.path.isfile(file):
                    os.remove(file)
                raise

            assert os.path.isfile(file)

            if os.path.isfile(file):
                with ZipFile(file) as z:
                    for png in z.namelist():
                        if os.path.isfile(png):
                            print(png)
                            os.remove(png)


def chop_main():
    def f_(file):
        with np.load(file) as npz:
            for f in npz.files:
                name = re.search(r'\d+', f)[0]
                try:
                    yield name, npz[f][:, :, :3].transpose(2, 0, 1) # (H, W, C) -> (C, H, W)
                except:
                    print(f)

    for file in tqdm(ut.globm('archive/*/result/out_*.npz')):
        out = file.replace('out', 'uvp')

        if os.path.isfile(out) and os.path.getsize(out) > 100 * 1048576:
            os.remove(file)

        else:
            try:
                np.savez_compressed(out, **dict(f_(file)))
                os.remove(file)

            except:
                print(file)
                if os.path.isfile(out):
                    os.remove(out)
                raise


    # with np.load(ut.globm('archive/*/result/uvp_*.npz')[0]) as npz:
    #     print(npz[npz.files[0]].shape)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comp', '-c', action='store_true', help='compress results')
    parser.add_argument('--zip', '-z', action='store_true', help='zip images')
    args = parser.parse_args()

    if args.comp:
        chop_main()

    if args.zip:
        with ut.chdir('archive'):
            zip_main()


if __name__ == '__main__':
    main()
