#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse

def visualize(df, rec_rates, snps, min, max, gargs, args):
    # df = df.sort_values(by=['stop', 'len'], ascending=False)
    df = df.sort_values(by=['start', 'len'], ascending=[True,False])

    y = [
        list(item) for item in zip(df['start'].to_list(), df['stop'].to_list())
    ]
    samples = df.id.values.tolist()
    samples = [x[4:] for x in samples]

    fig = plt.figure(figsize=(args.width, args.height))
    ax = fig.add_subplot(
        111,
        ylabel="sample",
        xlabel=gargs["xlabel"],
    )
    lines_color = '#e76f9b'
    # lines_color = '#8282e9'
    snp_color = '#FFF'
    text_color = '#000'
    rec_color = '#000'

    plt.setp(ax.get_xticklabels(), rotation=0, fontsize='small')
    ids = df.id.values.tolist()
    for idx in range(0, len(samples)):
        ax.plot(y[idx], [ids[idx], ids[idx]], c=lines_color, linewidth='2.5')

    ax.plot(rec_rates['pos'], rec_rates['rate'], c=rec_color, linewidth='1')

    plt.yticks(fontsize=3)
    plt.xlim(xmin=min, xmax=max)
    ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    for item in snps:
        ax.axvline(x=item[1] - 50,
                   ls='-',
                   c=snp_color,
                   linewidth='2')

        plt.text(item[1] + 500, len(y)/2, item[0], rotation=0, c=text_color, fontsize='x-large')

    ax.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.tick_params(axis='both', which='minor', labelsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=600)
    else:
        plt.show()

def load_df(name):
    df = pd.read_csv(name)

    sorted_start = df['start'].sort_values().tolist()
    sorted_stop = df['stop'].sort_values().tolist()

    # cutoff = int(len(sorted_start) * 0.2)
    cutoff = int(len(sorted_start) * 0.10)

    start_last = sorted_start[cutoff]
    stop_last = sorted_stop[-cutoff]

    df['start'] = df['start'].apply(lambda x: start_last if x < start_last else x)
    df['stop'] = df['stop'].apply(lambda x: stop_last if x > stop_last else x)


    df['len'] = df['stop'] - df['start']
    return df

def load_rec_rates(df, name, min, max):
    columns = ['chr','pos','rate','cm']
    rec_rates = pd.read_csv(name, header=None, names=columns, sep="\t")

    rec_rates = rec_rates[(rec_rates['pos'] >= min) & (rec_rates['pos'] <= max)]

    scale = df.shape[0] / (rec_rates['rate'].max() * 1.5)
    rec_rates['rate'] = rec_rates['rate'] * scale

    return rec_rates

pd.set_option('display.max_rows', 1000)

# C9ORF72
# snps = [("rs139185008", 27491944), ("rs117204439", 27607975), ("rs3849942", 27543283)]

# SOD1
# snps = [("SOD1*D91A", 31667290)]

# GSN

parser = argparse.ArgumentParser()
parser.add_argument('data', type=str)
parser.add_argument('-r', '--rec-rates', type=str, required=True)
parser.add_argument('-c', '--coords', type=str, required=True)
parser.add_argument('--gene', type=str, required=True)
parser.add_argument('--width', type=float, default=14.40)
parser.add_argument('--height', type=float, default=25.60)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df = load_df(args.data)

coords = args.coords
chr = coords.split(":")[0]
pos = coords.split(":")[1]
gene = args.gene

snps = [(gene, int(pos))]

min = df['start'].min() - 8000
max = df['stop'].max() + 10000

rec_rates = load_rec_rates(df, args.rec_rates, min, max)

graph_args = {
        "xlabel": chr,
        }

visualize(df, rec_rates, snps, min, max, graph_args, args)

