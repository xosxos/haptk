#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
import matplotlib

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

df = pd.read_csv(args.file)

df.sort_values(by=['allele'])
df = df[df["allele"] != 100]
df = df[df["markers"] > 10]
df = df[df["allele"] >= 10]

def rename_column(x):
    if x == 100:
        return "exp"
    else:
        return x

def select_color(alleles, x):
    set_obj = set(alleles)
    num = len(set_obj)
    allele = alleles[x]
    idx = list(set_obj).index(allele)
    palette = sns.color_palette("Set2", num)
    # print(allele, palette[idx])
    return palette[idx]

df_markers = df.copy()
# df_markers["allele"] = df_markers["allele"].apply(rename_column)

# fig, axs = plt.subplots(nrows=1,ncols=2)

# plt.bar(df_markers.allele, df_markers.markers, color ='maroon', width = 0.4)

fig = plt.figure(figsize=(28.80, 19.20))
ax = fig.add_subplot(
    111,
    # ylabel="sample",
    # xlabel=gargs["xlabel"],
)

markers = list(df_markers.markers)
alleles = list(df_markers.allele)
for idx in range(0, len(df_markers.allele)):
    color = select_color(alleles, idx)
    ax.plot([idx, idx], [0, markers[idx]], c = select_color(alleles, idx), linewidth='4')

ax.set_ylim([0, max(markers) + 100])
ax.set_xlim([-1, len(alleles) + 1])


def return_positions(alleles):
    counts = dict()
    for i in alleles:
      counts[i] = counts.get(i, 0) + 1

    count_list = list(map(lambda x: x[1], counts.items()))
    alleles = list(map(lambda x: x[0], counts.items()))
    positions = []

    for (i, count) in enumerate(count_list):
        positions.append((sum(count_list[0:i]) + count / 2) - 0.5)

    return (alleles[::2], positions[::2])


         

(labels, positions) = return_positions(alleles)
# positions = [(rect.get_x() + rect.get_width() / 2) for rect in ax.patches]
# labels = [round((rect.get_x() + rect.get_width() / 2) - 0.5) for rect in ax.patches]
ax.set_xticks(positions, labels)
ax.xaxis.set_tick_params(labelsize=11)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_ylabel('Segment length (markers)', fontsize=15)
 
# sns.barplot(y=df_markers.markers, x=df_markers.index, ax=axs[0]).set(title="")
# axs[0].set_xticklabels(df_markers['allele'])


# axs[0].tick_params(axis='both', which='major', labelsize=20)

# axs[0].spines['top'].set_visible(False)
# axs[0].spines['right'].set_visible(False)
# axs[0].spines['bottom'].set_visible(False)
# axs[0].spines['left'].set_visible(False)
# # axs[0].set_ylabel('Average segment length', fontsize=25)

# df = df[df["allele"] != 100]
# df = df[df["allele"] > 10]
# sns.histplot(data=df, x="allele", ax=axs[1], binwidth=1).set(title="")

# ax = sns.histplot(data=df, x="allele", binwidth=1)
# positions = [(rect.get_x() + rect.get_width() / 2) for rect in ax.patches]
# labels = [round((rect.get_x() + rect.get_width() / 2) - 0.5) for rect in ax.patches]
# axs[1].set_xticks(positions, labels)

# axs[1].tick_params(axis='both', which='major', labelsize=20)
# axs[1].spines['top'].set_visible(False)
# axs[1].spines['right'].set_visible(False)
# axs[1].spines['bottom'].set_visible(False)
# axs[1].spines['left'].set_visible(False)
# # axs[1].set_ylabel('MRCA estimate', fontsize=25)


plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
