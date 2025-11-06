#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker
import argparse

def visualize(fig, ax, df, rec_rates, snps, g, output, sample_colors, hap_start, hap_stop):
    # Use intra group sorting with the `group_order` key if several groups should be colored
    if sample_colors:
        df = df.sort_values(by=['group_order', 'len', 'start'], ascending=[False, False, True])
    else:
        # df = df.sort_values(by=['stop', 'len'], ascending=[False, False])
        # df = df.sort_values(by=['start', 'len'], ascending=[True,False])
        df = df.sort_values(by=['len', 'start'], ascending=[False, True])

    # Shared region for each sample in the format: [start, stop]
    shared_regions = [list(item) for item in zip(df['start'].to_list(), df['stop'].to_list())]

    # Sample names as UTF-8 strings
    samples = df.id.values.tolist()

    draw_horizontal_lines(ax, g, samples, shared_regions, hap_start, hap_stop)

    # Draw recombination rates
    if rec_rates is not None:
        ax.plot(rec_rates['pos'], rec_rates['rate'], c=g["rec_rates_color"], linewidth='1')
  
    draw_vertical_lines(plt, ax, g, snps, shared_regions)

    # Setup axis visibility and tick parameters etc
    setup_plot(ax, g)

    if output:
        plt.savefig(output, dpi=600)
    else:
        plt.show()


def draw_horizontal_lines(ax, g, samples, shared_regions, hap_start, hap_stop):
    line_color = g["shared_region_color"]
    hap_color = g["haplotype_color"]

   # Draw horizontal lines
    for idx in range(0, len(samples)):
        if sample_colors:
            lines_color = sample_colors[samples[idx]]

        x_start = shared_regions[idx][0]
        x_stop = shared_regions[idx][1]

        # If there is a haplotype, color it
        if hap_start is not None and hap_stop is not None:
            # Draw shared region
            ax.plot([x_start, x_stop], [samples[idx], samples[idx]], c=line_color, linewidth=g["line_width"])

            # Cut the haplotype if the shared region is smaller than the haplotype
            hap_start = max(x_start, hap_start)
            hap_stop = min(x_stop, hap_stop)

            # Draw haplotype region on top of the shared region
            ax.plot([hap_start, hap_stop], [samples[idx], samples[idx]], c=hap_color, linewidth=g["line_width"])
        else:
            ax.plot([x_start, x_stop], [samples[idx], samples[idx]], c=g["lines_color"], linewidth=g["line_width"])

    
def draw_vertical_lines(plt, ax, g, snps, shared_regions):
    count = 0

    for i, item in enumerate(snps):
        # Used to print label at 3 different heights to prevent overlapping
        if count >= 3:
            count = 1
        else:
            count += 1

        if i == 0:
            line_color = g["root_line_color"]
            text_color = g["root_label_color"]
            bg_color = g["root_bg_color"]
            height = len(shared_regions)/2
        else:
            line_color = g["snp_line_color"]
            text_color = g["snp_label_color"]
            bg_color = g["snp_bg_color"]

            # Change label height to avoid overlaps
            if count == 1:
                height = 13
            elif count == 2:
                height = 26
            else:
                height = 39
                
        # Draw line
        ax.axvline(x=item[1] - 50, ls='-', c=line_color, linewidth=g["snp_line_width"])

        # Draw associated text
        plt.text(
                 item[1] + g["snp_label_padding"],
                 height,
                 item[0],
                 rotation=0,
                 c=text_color,
                 fontsize=g["snp_font_size"],
                 backgroundcolor = bg_color
        )
    

def setup_plot(ax, g):
    if not g["show_ids"]:
        ax.get_yaxis().set_visible(False)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


    # Y axis params
    plt.yticks(fontsize=3)

    # ax.tick_params(axis='y', which='minor', labelsize=10)
    ax.tick_params(axis='y', which='major', labelsize=g["y_tick_size"])
    
    # X axis params
    print(g["xmin"] - g["x_axis_padding"])
    # plt.xlim(xmin=g["xmin"] - g["x_axis_padding"], xmax=g["xmax"] + g["x_axis_padding"])
    plt.xlim(g["xmin"] - g["x_axis_padding"], g["xmax"] + g["x_axis_padding"])
    ax.xaxis.label.set_size(15)

    # ax.tick_params(axis='x', which='minor', labelsize=10)
    ax.tick_params(axis='x', which='major', labelsize=g["x_tick_size"])

    # X tick formatter
    ax.get_xaxis().set_major_formatter( matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    plt.xlabel(g["xlabel"], fontdict={'size': g["x_label_size"]})
    plt.ylabel(g["ylabel"], fontdict={'size': g["y_label_size"]})

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


# ___MAIN CODE____
# 
pd.set_option('display.max_rows', 1000)

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
parser.add_argument('--line-width', type=int, default=5)
args = parser.parse_args()

df = load_df(args.data)

hap_start, hap_stop, rec_rates, gene = None, None, None, None

coords = args.coords
chr = coords.split(":")[0]
pos = coords.split(":")[1]

if args.gene:
    gene = args.gene
if args.hap_start:
    hap_start = args.hap_start
if args.hap_stop:
    hap_stop = args.hap_stop

xmin = df['start'].min() - 8000
xmax = df['stop'].max() + 10000

# Insert gene name as a special 'first' or 'root' snp to draw
if gene is not None:
    snps = [(gene, int(pos))]
else:
    snps = [("", int(pos))]

# Parse SNP lines (vertical lines) to draw 
if args.positions:
    with open(args.positions, 'r') as file:
        for line in file:
            line = line.strip()
            line = line.split(',')
            name = line[0]
            pos = int(line[1])
            if pos < xmax - 10000 and pos > xmin + 10000:
                snps.append((name, pos))

# Parse recombination rates
if args.rec_rates:
    rec_rates = load_rec_rates(df, args.rec_rates, xmin, xmax)


graph_args = {
        # General
        "xlabel": "",
        # "xlabel": chr,
        # "ylabel": "samples",
        "ylabel": "",
        "x_label_size": 10,
        "y_label_size": 10,
        "x_tick_size": 12,
        "y_tick_size": 10,
        "xmin": xmin,
        "xmax": xmax,
        "x_axis_padding": 300000,

        # Shared region lines
        "line_width": args.line_width,
        "haplotype_color": "#D0BDF0",
        "shared_region_color": '#df6a95',
        "group_colors": ["#866adb", "#f7708b", "#6ec9b8", "#ffe1bd", "#6ed4ff"],
        "show_ids": True,

        # Recombination rates
        "rec_rates_color": '#000',

        # Vertical haplotype lines
        "snp_label_padding": 110,
        "snp_line_width": 2,
        "snp_font_size": 'x-large',

        # The first item in the list can be colored differently
        "root_label_color": '#FFF',
        "root_line_color": '#FFF',
        "root_bg_color": '#808080',

        # The rest follow this color palette
        "snp_label_color": '#000',
        "snp_line_color": '#000',
        "snp_bg_color": '#FFF',
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

sample_colors = {}

for (vec, color) in zip(samples_to_tag, graph_args["group_colors"]):
    for sample in vec:
        sample_colors[sample] = color

# Add color information to the dataframe itself
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


# Initialize the figure
fig = plt.figure(figsize=(args.width, args.height))

ax = fig.add_subplot(
    111,
    # ylabel="sample",
    # xlabel=graph_args["xlabel"],
)

visualize(fig, ax, df, rec_rates, snps, graph_args, args.output, sample_colors, hap_start, hap_stop)

