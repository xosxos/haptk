#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse

def visualize(df, rec_rates, snps, min, max, gargs, args, sample_colors, hap_start, hap_stop):
    # df = df.sort_values(by=['stop', 'len'], ascending=[False, False])
    # df = df.sort_values(by=['start', 'len'], ascending=[True,False])
    df = df.sort_values(by=['len', 'start'], ascending=[False, True])

    if sample_colors:
        df = df.sort_values(by=['group_order', 'len', 'start'], ascending=[False, False, True])

    y = [list(item) for item in zip(df['start'].to_list(), df['stop'].to_list())]
    samples = df.id.values.tolist()

    fig = plt.figure(figsize=(args.width, args.height))

    ax = fig.add_subplot(
        111,
        ylabel="sample",
        xlabel=gargs["xlabel"],
    )

    hap_color = "#D0BDF0"
    lines_color = '#df6a95'
    snp_color = '#FFF'
    text_color = '#000'
    rec_color = '#000'

    ids = df.id.values.tolist()

    for idx in range(0, len(samples)):
        if sample_colors:
            lines_color = sample_colors[samples[idx]]

        x_start = y[idx][0]
        x_stop = y[idx][1]

        if hap_start is not None and hap_stop is not None:
            ax.plot([x_start, hap_start], [ids[idx], ids[idx]], c=lines_color, linewidth='5')
            ax.plot([hap_stop, x_stop], [ids[idx], ids[idx]], c=lines_color, linewidth='5')
            ax.plot([hap_start, hap_stop], [ids[idx], ids[idx]], c=hap_color, linewidth='5')
        else:
            ax.plot([x_start, x_stop], [ids[idx], ids[idx]], c=lines_color, linewidth='5')

    if rec_rates is not None:
        ax.plot(rec_rates['pos'], rec_rates['rate'], c=rec_color, linewidth='1')


    count = 0
    for i, item in enumerate(snps):

        padding = 110
        size = 'x-large'

        if count >= 3:
            count = 1
        else:
            count += 1

        if i == 0:
            color = "#000"
            ax.axvline(x=item[1] - 50, ls='-', c=color, linewidth='2')
            plt.text(item[1] + padding, len(y)/2, item[0], rotation=0, c="#fff", fontsize=size, backgroundcolor = "#000")
        else:
            ax.axvline(x=item[1] - 50, ls='-', c=snp_color, linewidth='2')

            if count == 1:
                height = 13
                plt.text(item[1] + padding, height, item[0], rotation=0, c=text_color, fontsize=size, backgroundcolor = "#fff")
            elif count == 2:
                height = 26
                plt.text(item[1] + padding, height, item[0], rotation=0, c=text_color, fontsize=size, backgroundcolor = "#fff")
            else:
                # height = len(y)/12
                height = 39
                # print(len(y))
                plt.text(item[1] + padding, height, item[0], rotation=0, c=text_color, fontsize=size, backgroundcolor = "#fff")
                
            

    setup_plot(ax)

    if args.output:
        plt.savefig(args.output, dpi=600)
    else:
        plt.show()

def setup_plot(ax):
    ax.get_yaxis().set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.yticks(fontsize=3)
    plt.xlim(xmin=min, xmax=max)
    ax.get_xaxis().set_major_formatter( matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    ax.tick_params(axis='y', which='minor', labelsize=10)
    ax.tick_params(axis='y', which='major', labelsize=10)

    ax.tick_params(axis='x', which='minor', labelsize=10)
    ax.tick_params(axis='x', which='major', labelsize=15)
    ax.xaxis.label.set_size(15)
    plt.tight_layout()

def load_df(name):
    df = pd.read_csv(name)

    # Cut x-axis at % of samples if so wanted
    percentage = 10

    if percentage and percentage != 0:
        sorted_start = df['start'].sort_values().tolist()
        sorted_stop = df['stop'].sort_values().tolist()
        cutoff = int(len(sorted_start) * (percentage / 100))
        start_last = sorted_start[cutoff]
        stop_last = sorted_stop[-cutoff]
        df['start'] = df['start'].apply(lambda x: start_last if x < start_last else x)
        df['stop'] = df['stop'].apply(lambda x: stop_last if x > stop_last else x)


    # Add length column
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
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('-r', '--rec-rates', type=str)
parser.add_argument('-c', '--coords', type=str, required=True)
parser.add_argument('--gene', type=str)
parser.add_argument('-x', '--width', type=float, default=14.40)
parser.add_argument('-y', '--height', type=float, default=25.60)
parser.add_argument('--hap-start', type=int, required=False)
parser.add_argument('--hap-stop', type=int, required=False)
parser.add_argument('-o', '--output', type=str, required=False)
parser.add_argument('--positions', type=str, required=False)
args = parser.parse_args()

df = load_df(args.data)

hap_start = None
hap_stop = None
rec_rates = None
gene = None

coords = args.coords
chr = coords.split(":")[0]
pos = coords.split(":")[1]

if args.gene:
    gene = args.gene
if args.hap_start:
    hap_start = args.hap_start
if args.hap_stop:
    hap_stop = args.hap_stop

min = df['start'].min() - 8000
max = df['stop'].max() + 10000

if gene is not None:
    snps = [(gene, int(pos))]
else:
    snps = [("", int(pos))]
    
if args.positions:
    with open(args.positions, 'r') as file:
        for line in file:
            line = line.strip()
            line = line.split(',')
            name = line[0]
            pos = int(line[1])
            if pos < max - 10000 and pos > min + 10000:
                snps.append((name, pos))

if args.rec_rates:
    rec_rates = load_rec_rates(df, args.rec_rates, min, max)

graph_args = {
        "xlabel": chr,
        }


# Create list of samples to tag
samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

# kyme, kanta-hame, karjala, paijat-hame, pirkanmaa
colors=["#866adb", "#f7708b", "#6ec9b8", "#ffe1bd", "#6ed4ff"]
# colors=["#c34757", "#9a69b1", "#4071ab", "#90a720", "#d75bae"]

sample_colors = {}
for (vec, color) in zip(samples_to_tag, colors):
    for sample in vec:
        sample_colors[sample] = color

if sample_colors:
    df['color'] = df.apply(lambda row: sample_colors[row.id], axis = 1)

    tuples = []
    for color in df.color.unique():
        color_df = df.loc[df['color'] == color]
        sum = color_df.len.sum()
        tuples.append((color, sum/color_df.shape[0]))

    tuples = sorted(tuples, key=lambda x: float(x[1]))

    order = {}
    for i, elem in enumerate(tuples):
        order[elem[0]] = i

    # Lowest order number has lowest haplotype sharing on avg
    df['group_order'] = df.apply(lambda row: order[row.color], axis = 1)

visualize(df, rec_rates, snps, min, max, graph_args, args, sample_colors, hap_start, hap_stop)

