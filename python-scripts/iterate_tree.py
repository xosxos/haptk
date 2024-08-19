import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter

# Import the HAPTK python library
import haptk

### Script logic
parser = argparse.ArgumentParser()
parser.add_argument('hst', type=str)
parser.add_argument('--df', type=str)
parser.add_argument('--min-size', type=int, default=1)    
parser.add_argument('--hard-cut', action="store_true")    
parser.add_argument('--ids', nargs="+", type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

# Read an .hst.gz file
hst = haptk.read_hst(args.hst)

df = pd.read_csv(args.df, sep=",")

# print(hst.metadata)
samples = hst.metadata["samples"]

mask = df['id'].isin(samples)
df = df[mask]


calculate_dur = True
# Calculate disease duration for those not deceased
if calculate_dur:
    # Up to date survival data is not available on 83 samples
    df["status"] = df.apply(lambda row: 0 if np.isnan(row["dur"]) else 1, axis=1)

    # Calculate duration from the year of onset (data has been last checked in 2021, thus subtract from 2021)
    check_yoo = lambda row: np.nan if np.isnan(row['yoo']) else 2021 - row['yoo']

    df["dur"] = df.apply(lambda row: check_yoo(row) if np.isnan(row["dur"]) else row["dur"], axis=1)
    df = df.drop("yoo", axis = 1)
    df.dropna(inplace=True)
else:
    df["status"] = 1
    df = df.drop("yoo", axis = 1)
    df.dropna(inplace=True)


df["site"] = df["site"].apply(lambda x: 0 if x == "SPINAL" else 1)
    
# Code sex as male=0 and female=1
df["sex"] = df["sex"].apply(lambda x: 0 if x == "M" else 1)

# print(df)

# df = df.drop("id", axis = 1)

# # Fit the model
# cph = CoxPHFitter()

# cph.fit(df, duration_col = 'dur', event_col = 'status')

# len = len(cph.summary['coef'])

# print(f"covariate,coef,exp_coef,se_coef,coef_lower_95,coef_upper_95,exp_coef_lower_95,exp_coef_upper_95,cmp_to,z,p,minus_log2_p")
# for i in range(0,len):
#     covariate = cph.summary['coef'].index[i]
#     if covariate == "aoo":
#         coef = round(cph.summary['coef'].iloc[i], 6)
#         exp_coef = round(cph.summary['exp(coef)'].iloc[i], 6)
#         se_coef = round(cph.summary['se(coef)'].iloc[i], 6)
#         coef_lower_95 = round(cph.summary['coef lower 95%'].iloc[i], 6)
#         coef_upper_95 = round(cph.summary['coef upper 95%'].iloc[i], 6)
#         exp_coef_lower_95 = round(cph.summary['exp(coef) lower 95%'].iloc[i], 6)
#         exp_coef_upper_95 = round(cph.summary['exp(coef) upper 95%'].iloc[i], 6)
#         cmp_to = round(cph.summary['cmp to'].iloc[i], 6)
#         z = round(cph.summary['z'].iloc[i], 6)
#         p = cph.summary['p'].iloc[i]
#         p = f"{p:.4g}"
#         minus_log2_p = round(cph.summary['-log2(p)'].iloc[i], 6)

#         print(f"{covariate},{coef},{exp_coef},{se_coef},{coef_lower_95},{coef_upper_95},{exp_coef_lower_95},{exp_coef_upper_95},{cmp_to},{z},{p},{minus_log2_p}")

# print(f"covariate,coef,exp_coef,se_coef,coef_lower_95,coef_upper_95,exp_coef_lower_95,exp_coef_upper_95,cmp_to,z,p,minus_log2_p")
# for i in range(0,len):
#     covariate = cph.summary['coef'].index[i]
#     coef = round(cph.summary['coef'].iloc[i], 6)
#     exp_coef = round(cph.summary['exp(coef)'].iloc[i], 6)
#     se_coef = round(cph.summary['se(coef)'].iloc[i], 6)
#     coef_lower_95 = round(cph.summary['coef lower 95%'].iloc[i], 6)
#     coef_upper_95 = round(cph.summary['coef upper 95%'].iloc[i], 6)
#     exp_coef_lower_95 = round(cph.summary['exp(coef) lower 95%'].iloc[i], 6)
#     exp_coef_upper_95 = round(cph.summary['exp(coef) upper 95%'].iloc[i], 6)
#     cmp_to = round(cph.summary['cmp to'].iloc[i], 6)
#     z = round(cph.summary['z'].iloc[i], 6)
#     p = cph.summary['p'].iloc[i]
#     p = f"{p:.4g}"
#     minus_log2_p = round(cph.summary['-log2(p)'].iloc[i], 6)

#     print(f"{covariate},{coef},{exp_coef},{se_coef},{coef_lower_95},{coef_upper_95},{exp_coef_lower_95},{exp_coef_upper_95},{cmp_to},{z},{p},{minus_log2_p}")

# Create list of samples to tag
samples_to_tag = []
if args.ids:
    for file in args.ids:
        ids = []
        file = open(file, 'r')
        for line in file.readlines():
            ids.append(line.strip())

        samples_to_tag.append(ids)

def optimizer(hst, node_name, indexes, df):
    # print(node_name)
    if node_name != "0":
        if hst.metadata['ploidy'] == "Diploid":
            samples = hst.metadata["samples"]
            samples = samples + hst.metadata["samples"]
        else:
            samples = hst.metadata["samples"]
        
        names = []
        for i in indexes:
            name = samples[i]
            names.append(name)

        df["gt"] = df["id"].apply(lambda id: names.count(id))

        df = df.drop("id", axis = 1)


        # print(list(df["gt"]))
        hets = list(df["gt"]).count(1)
        alt_homs = list(df["gt"]).count(2)
        if hets + alt_homs < 5:
            return 1.0

        # Fit the model
        cph = CoxPHFitter()

        cph.fit(df, duration_col = 'dur', event_col = 'status')

        coef = round(cph.summary['coef'].iloc[5], 6)
        p = cph.summary['p'].iloc[5]
        # p = f"{p:.4g}"
        print(node_name, len(indexes), coef, p)

        return float(p)
    else:
        return 1.0



# # Render the tree
hst.iterate_tree(df, optimizer, args.output, w=4600, h=4600, to_tag=samples_to_tag, min_size=args.min_size, hard_cut=args.hard_cut)




