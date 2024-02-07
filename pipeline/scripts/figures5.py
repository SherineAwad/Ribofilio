#! /usr/bin/env python
import sys
import argparse
import screed
import math
import numpy as np
import matplotlib.pyplot as plt
import csv





def main():
   fig = plt.figure(figsize =(15, 15))
   ax = fig.add_subplot(111)
   rbs = np.array([0.0051,0.0102,0.0100,0.0073])
   rcs = np.array([0.0003, 0.0006, 0.0006, 0.0004])
   ds = np.array(['D1', 'D2', 'D3', 'D4'])


   D=['D1','D2','D3','D4']
   r=['$r_b$','$r_c$']
   pos = np.arange(len(D))
   rbs = np.array([0.0051,0.0102,0.0100,0.0073])
   rcs = np.array([0.0003, 0.0006, 0.0006, 0.0004])

   plt.bar(pos,rbs,color='blue',edgecolor='black')
   plt.bar(pos,rcs,color='grey',edgecolor='black')
   plt.xticks(pos, D, fontsize =20)
   plt.yticks(fontsize=20)
   plt.xlabel('Datasets', fontsize=20)
   plt.ylabel('$r_b$ and  $r_c$', fontsize=20)
   plt.title('Datasets dropoff rate',fontsize=20)
   plt.legend(r,loc=2, fontsize=20)
   plt.savefig("datasets.png", format="png")


if __name__ == '__main__':
    main()

