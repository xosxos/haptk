#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


parser = argparse.ArgumentParser()
parser.add_argument('data', type=str)
parser.add_argument('--ctrls', type=str)
parser.add_argument('--window-length', type=int, default=101)
parser.add_argument('--polyorder', type=int, default=3)
parser.add_argument('--ymax', type=int, default=4)
parser.add_argument('--width', type=int, default=2560)
parser.add_argument('--height', type=int, default=1440)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df: pd.DataFrame = pd.read_csv(args.data, usecols=['POS', 'mrca'])

df["mrca"] = np.log10(df["mrca"])

mrca2 = savgol_filter(df.mrca, args.window_length, args.polyorder)

fig = plt.figure(figsize=(args.width / 100, args.height / 100))
ax = fig.add_subplot(111)

# ax.plot(df.pos, df.mrca, linestyle='-', color="#ec66e6", linewidth=4, alpha=.25)
ax.plot(df["POS"], df["mrca"], linestyle='-', color="#f62323", linewidth=3, alpha=.5)

# std = np.std(df.mrca)
# mean = np.mean(df.mrca)
# sig_line = mean - 4 * std
# plt.axhline(y = sig_line, color="#222", linestyle=(0, (3, 5, 1, 5)), linewidth=4, alpha=.5)

if args.ctrls:
    ctrl_df: pd.DataFrame = pd.read_csv(args.ctrls, usecols=['POS', 'mrca'])
    ctrl_mrca2 = savgol_filter(ctrl_df.mrca, args.window_length, args.polyorder)
    ax.plot(ctrl_df["POS"], ctrl_mrca2, color='#484848', linewidth=4, alpha=1)

ax.plot(df["POS"], mrca2, color='black', linewidth=3.5, alpha=1)

ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)

ax.get_xaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.xlim(xmin=df["POS"].min(), xmax=df["POS"].max())
plt.ylim(ymin=1.0, ymax=args.ymax)


plt.title(args.data)

plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()


