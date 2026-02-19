#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
import matplotlib.ticker as ticker
import matplotlib as mpl

from matplotlib.ticker import MultipleLocator
# various formatting parameters
label_fontsize = 10
tick_fontsize = 10
linewidth = 1
major_xtick_length = 15
minor_xtick_length = 7
major_ytick_length = 7
minor_ytick_length = 0

mpl.rcParams['figure.dpi'] = 400
mpl.rcParams['font.weight'] = 'normal'
mpl.rcParams['axes.linewidth'] = linewidth
mpl.rcParams['lines.linewidth'] = linewidth
mpl.rcParams['xtick.labelsize'] = tick_fontsize
mpl.rcParams['ytick.labelsize'] = tick_fontsize
mpl.rcParams['xtick.major.width'] = linewidth
mpl.rcParams['ytick.major.width'] = linewidth
mpl.rcParams['xtick.minor.width'] = linewidth
mpl.rcParams['ytick.minor.width'] = linewidth

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df = pd.read_csv(args.file)

df.sort_values(by=['allele'])
df = df[df["allele"] != 100]
# df = df[df["markers"] > 10]
# df = df[df["allele"] >= 10]

def rename_column(x):
    if x == 100:
        return "exp"
    else:
        return x

def select_color(alleles, idx):
    set_obj = set(alleles)
    num = len(set_obj)
    allele = alleles[idx]
    idx = list(set_obj).index(allele)
    palette = sns.color_palette("Set2", num)
    # print(allele, palette[idx])
    return palette[idx]

df_markers = df.copy()
# df_markers["allele"] = df_markers["allele"].apply(rename_column)

# fig, axs = plt.subplots(nrows=1, ncols=2)

# plt.bar(df_markers.allele, df_markers.markers, color ='maroon', width = 0.4)

fig = plt.figure(figsize=(12.0, 4.0))
ax = fig.add_subplot(
    111,
    # ylabel="sample",
    # xlabel=gargs["xlabel"],
)

markers = list(df_markers.markers)
alleles = list(df_markers.allele)
sex = list(df_markers.parent_sex)
for idx in range(0, len(df_markers.allele)):
    # color = select_color(alleles, idx)
    if sex[idx] == "F":
        # color = "#FFCF42"
        color = "#B32656"
    else:
        color = "#0073B2"

    # Make the alleles with length=0 also visible
    if markers[idx] == 0:
        y_range = [-0.15, 0.15]
    else:
        if markers[idx] < 0:
            y_range = [0.15, markers[idx]]
        else:
            y_range = [-0.15, markers[idx]]
            

            
        
    ax.plot([idx, idx], y_range, c = color, linewidth='5', solid_capstyle='butt')

# ax.set_ylim([min(markers), max(markers)])


def return_positions(alleles):
    counts = dict()

    # get allele counts
    for i in alleles:
      counts[i] = counts.get(i, 0) + 1

    count_list = list(map(lambda x: x[1], counts.items()))
    alleles = list(map(lambda x: x[0], counts.items()))
    major_positions = []
    minor_positions = []

    for (i, count) in enumerate(count_list):
        minor_positions.append((sum(count_list[0:i]) + count / 2) - 0.5)
        major_positions.append((sum(count_list[0:i])) - 0.5)

    # return (alleles[::2], positions[::2])
    return (alleles, major_positions, minor_positions)



# CONFIGURE X TICKS
(labels, major_positions, minor_positions) = return_positions(alleles)

ax.set_xticks(major_positions, [], minor=False)
ax.set_xticks(minor_positions, labels, minor=True)

# Get current ticks
# Adjust every other tick label
for i, tick in enumerate(ax.xaxis.get_minor_ticks()):
    # Hide every other tick on the top and bottom
    if i % 2 == 0:  # For every other tick
        tick.label2.set(alpha=0, text="", label="", color = "white")

        # plt.gca().get_xticklabels()[i].set_verticalalignment('bottom')
    else:
        tick.label1.set(alpha=0, text="", label="", color = "white")
        # plt.gca().get_xticklabels()[i].set_verticalalignment('top')


# CONFIGURE Y TICKS
# plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MultipleLocator(2))

ax.set_ylim([-4, 8])
ax.set_xlim([-1, len(alleles) + 1])

ax.set_ylabel('C H A N G E   ( R E P E A T S )', fontsize=label_fontsize)
ax.set_xlabel('P A R E N T   A L L E L E', fontsize=label_fontsize, labelpad = 10)

ax.tick_params('x', which='both', bottom=True, top=True, labeltop=True, direction='in', labelsize=tick_fontsize)
ax.tick_params(axis='x', which='minor', tick1On=False, tick2On=False)
ax.tick_params('y', left=True, right=True, direction='in', labelright=True, labelsize=tick_fontsize)

# plt.grid(True, color="#252525", linestyle=(0, (1, 10)), axis='x')
# plt.grid(True, color="#252525", linestyle=(0, (10, 10)), axis='x')
plt.grid(True, color="#555555", axis='y')
# plt.grid(True, color="#555555", axis='x')
# plt.grid(linestyle="--")

xticks, _ = plt.xticks()
for x0, x1 in zip(xticks[::2], xticks[1::2]):
    plt.axvspan(x0, x1, color='black', alpha=0.05, zorder=0)
plt.xticks(xticks)  # force the same yticks again

plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
