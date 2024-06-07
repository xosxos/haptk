#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data', type=str)
parser.add_argument('--height', type=int, default=1440)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df = pd.read_csv(args.data)

# print(df)

fig, axs = plt.subplots(nrows=1,ncols=2)
# fig, ax = plt.subplots(nrows=1,ncols=1)

p = sns.color_palette("Set2")
cols=[p[0], p[5]]
cols=["#74e3e3", "#dc6ca6"]

df["Length (Mb)"] = df.length / 1000000
print(df["Length (Mb)"].describe())
print(df["Length (Mb)"].median())

sns.histplot(y=df["Length (Mb)"], ax=axs[0], bins=90, element="step", fill=True, color="#5f82b5").set(title="")
# sns.histplot(x=df.length, ax=ax, kde=True, bins=70).set(title="")
# sns.violinplot(y=df.length, ax=axs[1], cut=0).set(title="")
sns.violinplot(data=df, y=df["Length (Mb)"], cut=0, ax=axs[1], saturation=0.80, inner_kws=dict(box_width=10, whis_width=5, color="black"))
# sns.boxplot(y=df.length, ax=axs[2]).set(title="")
# sns.histplot(x=df.markers, ax=axs[1], kde=True).set(title="")

# axs[0].set_ylabel('Avg. segment length (markers)', fontsize=20)
# axs[1].set_ylabel('MRCA estimate (generations)', fontsize=20)

# axs[0].tick_params(axis='both', which='major', labelsize=20)
# axs[1].tick_params(axis='both', which='major', labelsize=20)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_visible(False)

# axs[0].get_yaxis().set_visible(False)
# axs[1].get_yaxis().set_visible(False)
axs[0].spines['top'].set_visible(False)
axs[1].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[1].spines['right'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[1].spines['left'].set_visible(False)
# axs[0].get_legend().remove()
# axs[1].get_legend().remove()

# axs[0].set_xlabel('Group', fontsize=25)
# axs[1].set_xlabel('Group', fontsize=25)


# plt.tight_layout()

plt.ticklabel_format(style='plain', axis='y')

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
