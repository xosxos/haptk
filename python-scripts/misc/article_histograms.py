#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--affy', type=str)
parser.add_argument('--illu', type=str)
parser.add_argument('--height', type=int, default=1440)
parser.add_argument('-o', '--output', type=str, required=False)
args = parser.parse_args()

illu = pd.read_csv(args.illu)
affy = pd.read_csv(args.affy)

illu['cohort'] = 'Tampere ALS-FTD and CTRL cohorts'
affy['cohort'] = 'Helsinki ALS'

df = pd.concat([illu, affy], ignore_index = True)
print(df)

fig, axs = plt.subplots(nrows=1,ncols=2)

p = sns.color_palette("Set2")
cols=[p[0], p[5]]
cols=["#74e3e3", "#dc6ca6"]

sns.barplot(y=df.avg_markers, x=df.group, hue=df.cohort, ax=axs[0], palette=cols).set(title="")
sns.barplot(y=df.mrca, x=df.group, hue=df.cohort, ax=axs[1], palette=cols).set(title="")

axs[0].set_ylabel('Avg. segment length (markers)', fontsize=20)
axs[1].set_ylabel('MRCA estimate (generations)', fontsize=20)

axs[0].tick_params(axis='both', which='major', labelsize=20)
axs[1].tick_params(axis='both', which='major', labelsize=20)

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
axs[0].get_legend().remove()
axs[1].get_legend().remove()

# axs[0].set_xlabel('Group', fontsize=25)
# axs[1].set_xlabel('Group', fontsize=25)


# plt.tight_layout()

if args.output:
    plt.savefig(args.output)
else:
    plt.show()
