# -*- coding: utf-8 -*-
import pandas as pd
import sys

args = sys.argv

data = pd.DataFrame()

with open(args[1]) as line:
    for i in line:
        list1 = i.strip().split("\t")
        file_n = str(list1[0]) + ".starded.csv"
        data1 = pd.read_csv(file_n)
        data1 = data1[(data1['MutationFrequency'] < 0.9) & (data1['MutationFrequency'] > 0.1)]
        if data.empty:
            data = data1
        else:
            data = pd.concat([data, data1], ignore_index=True)

data.to_csv('GATK_information.csv', index=False)
