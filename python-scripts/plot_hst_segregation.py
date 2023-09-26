#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import scipy.stats as stats

parser = argparse.ArgumentParser()
parser.add_argument('data', type=str)
parser.add_argument('--window-length', type=int, default=201)
parser.add_argument('--polyorder', type=int, default=3)
parser.add_argument('--width', type=int, default=2560)
parser.add_argument('--height', type=int, default=1440)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df: pd.DataFrame = pd.read_csv(args.data)

fig = plt.figure(figsize=(args.width / 100, args.height / 100))
ax = fig.subplots(2, 1)

savgol_sum = savgol_filter(df.bp_len, args.window_length, args.polyorder)
ax[0].plot(df["POS"], df.bp_len, linestyle='-', color="blue", linewidth=3, alpha=.20)
ax[0].plot(df["POS"], savgol_sum, linestyle='-', color="blue", linewidth=3, alpha=.90)

# df['freq'] = df['ctrls_freq'].apply(lambda x: 4 if x == 0 else -np.log10(x))
# ax[2].scatter(df.pos, df.freq, color="orange", s=2, alpha=.90)

ax[1].plot(df["POS"], df["n"], color="red", alpha=.90)


# std = np.std(df.value)
# mean = np.mean(df.value)
# sig_line = mean + 4 * std
# plt.axhline(y = sig_line, color="#222", linestyle=(0, (3, 5, 1, 5)), linewidth=4, alpha=.5)

# ax.plot(df.pos, savgol, color='black', linewidth=3, alpha=1)

for i, axis in enumerate(ax):
    ax[i].set_xlim(xmin=df["POS"].min(), xmax=df["POS"].max())
    ax[i].spines.right.set_visible(False)
    ax[i].spines.top.set_visible(False)
    ax[i].spines.bottom.set_visible(False)

    ax[i].get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    # ax[i].axvline(x=27574000, ls='--', c='#000', label="c9orf72", linewidth='1')
    # ax[i].axvline(x=121302000, ls='--', c='#000', label="gsn", linewidth='1')

    if i == 0:
        ymax = max(10, df.bp_len.max())
        ax[i].set_ylim(ymin=0, ymax=ymax)
        ax[i].title.set_text(args.data)

    if i == 1:
        ymax = max(10, df["n"].max())

        try:
            ax[i].set_ylim(ymin=0, ymax=ymax)
        except:
            ax[i].set_ylim(ymin=0, ymax=30)
        # ax[i].get_xaxis().set_visible(False)
        
plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()


