
import pandas as pd

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from lifelines import KaplanMeierFitter
import argparse

# pd.set_option("display.max_rows", None, "display.max_columns", None)

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str)
parser.add_argument('ids', type=str)
parser.add_argument('pheno', type=str)
args = parser.parse_args()

df = pd.read_csv(args.file)


lines = []
with open(args.ids) as file:
    lines = [line.rstrip() for line in file]

df["maj"] = df["id"].apply(lambda x: 0 if x in lines else 1)
df["site"] = df["site"].apply(lambda x: "LE" if x == "LIMB" else x)
df["site"] = df["site"].apply(lambda x: "BU" if x == "RESP" else x)

dummies_site = pd.get_dummies(df["site"], prefix = 'site')
df = pd.concat([df, dummies_site], axis = 1)

if args.pheno == "dur":
    df = df.dropna()

# x = c9df[['gender', 'site', 'gt']]
x = df[['maj']]
x = pd.get_dummies(data=x, drop_first=True)
y = df[args.pheno]
x = sm.add_constant(x)  # adding a constant
model = sm.OLS(y, x).fit()
print(model.pvalues[1])

# fig, axs = plt.subplots(3, 2)

# sns.regplot(data=df, y=args.pheno, x="maj", x_jitter=.1, ax=axs[0, 0])
# sns.boxplot(data=df, y=args.pheno, x="maj", ax=axs[0, 1])

# fig.tight_layout()
# plt.savefig(f"{args.pheno}_c9.png")

crosstab = pd.crosstab(df["maj"],  df["site_BU"], normalize='index')
print(crosstab)

crosstab = pd.crosstab(index=df["maj"],  columns=df["site_BU"])
results = stats.fisher_exact(crosstab)
print(results[1])
# crosstab = pd.crosstab(index=df["maj"],  columns=df["site_UE"])
# results = stats.fisher_exact(crosstab)
# print(results[1])
# crosstab = pd.crosstab(index=df["maj"],  columns=df["site_LE"])
# results = stats.fisher_exact(crosstab)
# print(results[1])


