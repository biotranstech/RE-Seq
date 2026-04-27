# -*- coding: utf-8 -*-
import pandas as pd
import sys


args = sys.argv

sample1 = args[1]
sample2 = args[2]
sample3 = args[3]

sample1 = pd.read_csv(sample1)
sample2 = pd.read_csv(sample2)
sample2 = sample2.iloc[:,5:]

result = pd.concat([sample1, sample2], axis=1)

result.to_csv(sample3 + '.starded.csv', index=False)


