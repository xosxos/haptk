import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from scipy.signal import savgol_filter
import re
import seaborn as sns

def draw_lines(df, angle, cum_angle, color):
    # Set the coordinates limits
    lowerLimit = 0

    # Compute max and min in the dataset
    max = df['savgol'].max()

    # Let's compute heights: they are a conversion of each item value in those new coordinates
    # In our example, 0 in the dataset will be converted to the lowerLimit (10)
    # The maximum will be converted to the upperLimit (100)
    slope = (max - lowerLimit) / max
    heights = slope * df["savgol"] + lowerLimit

    # Compute the width of each bar. In total we have 2*Pi = 360°
    width = angle / len(df.index)

    # Compute the angle each bar is centered on:
    indexes = list(range(1, len(df.index)+1))
    angles = [cum_angle + element * width for element in indexes]

    # Draw Savitzky line
    bars = ax.plot(
        angles, 
        heights,
        linewidth=6,
        color=color
    )

    print(max, lowerLimit)
    slope = (max - lowerLimit) / max
    heights = slope * df[column] + lowerLimit

    # Compute the width of each bar. In total we have 2*Pi = 360°
    width = angle / len(df.index)

    # Compute the angle each bar is centered on:
    indexes = list(range(1, len(df.index)+1))
    angles = [cum_angle + element * width for element in indexes]

    # Draw raw data
    bars = ax.plot(
        angles, 
        heights,
        linewidth=1,
        alpha = 0.2,
        color=color
    )


def deserved_radians(df_len, total_len):
    return (2 * np.pi) * (df_len / total_len)

def centromeres_hg38(chr):
    chrs = {
        "chr1": (121700000, 125100000),
        "chr2": (91800000, 96000000),
        "chr3": (87800000, 94000000),
        "chr4": (48200000, 51800000),
        "chr5": (46100000, 51400000),
        "chr6": (58500000, 62600000),
        "chr7": (58100000, 62100000),
        "chr8": (43200000, 47200000),
        "chr9": (42200000, 45500000),
        "chr10": (38000000, 41600000),
        "chr11": (51000000, 55800000),
        "chr12": (33200000, 37800000),
        "chr13": (16500000, 18900000),
        "chr14": (16100000, 18200000),
        "chr15": (17500000, 20500000),
        "chr16": (35300000, 38400000),
        "chr17": (22700000, 27400000),
        "chr18": (15400000, 21500000),
        "chr19": (24200000, 28100000),
        "chr20": (25700000, 30400000),
        "chr21": (10900000, 13000000),
        "chr22": (13700000, 17400000),
        "chrX": (58100000, 63800000),
        "chrY": (10300000, 10600000),
    }
    return chrs[chr]

parser = argparse.ArgumentParser()
parser.add_argument('-d','--dfs', nargs='+', help='<Required> Set flag', required=True)
parser.add_argument('-c','--column', type=str, help='<Required> Set flag', required=True)
parser.add_argument('--window-length', type=int, default=400)
parser.add_argument('--polyorder', type=int, default=1)
parser.add_argument('--ymax', type=int, default=4)
parser.add_argument('--width', type=int, default=2560)
parser.add_argument('--height', type=int, default=1440)
parser.add_argument('-o', '--output', type=str, required=False)

args = parser.parse_args()

column = args.column
print(args.dfs)
dfs = []
chrs = []

plt.figure(figsize=(20,20))
ax = plt.subplot(111, polar=True)

for filename in args.dfs:
    x = re.search("chr.*?_", filename).group(0)
    chr = x[:-1]

    df = pd.read_csv(filename, usecols=['pos', column, 'centromere'])

    df = df[df.centromere != True]

    # Drop markers from centromeres to remove centromeric effects
    df = df[~df["pos"].between(centromeres_hg38(chr)[0]-500_000, centromeres_hg38(chr)[1]+500_000)]
    print(len(df))

    if chr == "chr13":
        print(df)

    if column == 'mrca':
        df[column] = np.log10(df[column])
        plt.ylim(ymin=1.0, ymax=4)

    df[column] = df[column]
    df["savgol"] = savgol_filter(df[column], args.window_length, args.polyorder)
    dfs.append(df)
    chrs.append(chr)

# plt.axis('off')

colors = sns.color_palette("colorblind", 28)
colors = colors[4:28]

total_len = sum([len(df) for df in dfs])

cum_radians = 0;

x_tick_positions = []

for (df, color) in zip(dfs, colors):
    df_radians = deserved_radians(len(df), total_len)
    draw_lines(df, df_radians, cum_radians, color)
    x_tick_positions.append(cum_radians + df_radians / 2)
    cum_radians += df_radians


ax.set_xticks(x_tick_positions)
ax.set_xticklabels(chrs)
ax.tick_params(axis='both', which='major', pad=16)  # move the tick labels
ax.grid("lightgray", linestyle='--', linewidth=0.5)

ax.spines['polar'].set_visible(False)

plt.gca().invert_yaxis()
plt.setp(ax.get_xticklabels(), rotation=0, fontsize='28')
plt.setp(ax.get_yticklabels(), rotation=0, fontsize='25')
plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
