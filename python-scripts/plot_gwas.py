#!/bin/python3

import gwaslab as gl
import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data', type=str)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

df: pd.DataFrame = pd.read_csv(args.data)
print(df)

stats = gl.Sumstats(df, fmt="saige", snpid="MarkerID")
stats.basic_check()
stats.plot_mqq(anno=True, save=args.output, sig_level=1e-8,)

