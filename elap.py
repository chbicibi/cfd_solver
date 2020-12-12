from myutils import *
import myutils as ut
ut._import('numpy', 'np', table=globals())
ut._import('pandas', 'pd', table=globals())
ut._import('matplotlib', table=globals())
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter


def main1():
    odir = 'test'
    with ut.chdir(odir):
        times = sorted((os.path.getmtime(file), file)
                       for file in ut.iglobm('__raw__/*.plt')
                       if int(re.search(r'\d+', file)[0]) < 1000)
        a = np.array([t2[0] - t1[0] for t1, t2 in zip(times[:-1], times[1:])])
        list(map(print, (t[1] for t in times)))
        plt.plot(a)
        plt.xlim(0, len(a)-1)
        # plt.ylim(0, 1000)
        plt.show()


def main():
    plt.figure()
    for file in ut.iglobm('*{ada*,o1.7-pblc*}/nitr.csv'):
        if 'fix' in file or 'colab' in file:
            continue
        label = os.path.dirname(file)

        df = pd.read_csv(file)
        a = np.asarray(df, dtype=np.float64)

        num = 1000
        a[:, 1] = np.asarray([a[max(i-num+1, 0):i+1, 1].mean()
                             for i in range(len(a))])

        plt.plot(*a.T, lw=1, label=label)

    plt.xlim(0, 20000)
    plt.ylim(0, 100)
    # plt.yscale('log')
    plt.legend()
    plt.xlabel('cycle')
    plt.ylabel('itr')
    # plt.show()
    png = 'nitr-ada.png'
    plt.savefig(png, bbox_inches='tight', pad_inches=0.1, dpi=600)


if __name__ == '__main__':
    main()
