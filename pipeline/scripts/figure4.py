#! /usr/bin/env python
import argparse
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt


def get_R (DS):
    df = []
    normalised_dropoff = [] 
    for line in open(DS):
          df.append(line.strip()) 
    dropoff = np.asarray(df, dtype = float)
    mean = np.mean(dropoff)
    std = np.std(dropoff)
    #print(dropoff)
    for i in range(0, len(dropoff)) :
         normalised_dropoff.append( ( abs(dropoff[i]) - mean )/std ) 
         #print(normalised_dropoff[i]) 
    return normalised_dropoff
def main():
    n1 = get_R("d1_r.txt")
    n2 = get_R("d2_r.txt")
    n3 = get_R("d3_r.txt") 
    n4 = get_R("d4_r.txt")
    n5 = get_R("d5_r.txt") 
    n6 = get_R("d6_r.txt") 
    n7 = get_R("d7_r.txt") 
    n8 = get_R("d8_r.txt") 
    gene_length = [0,1,2,3,4,5]
    print(len(n1), len(n2), len(n3), len(n4), len(n5), len(n6), len(n7), len(n8)) 
    
    fig = plt.figure(figsize=(24,24))
    plt.ylabel("Normalised Drop-off rate", fontsize=22) 
    plt.xlabel("Gene Length interval", fontsize=22) 
    labels = ["]500-1000]","]1000-2000]","]2000-3000]","]3000-4000]","]4000-5000]" ,">5000"]
    plt.xticks(gene_length, labels, rotation='vertical', fontsize=22)
    plt.yticks(fontsize=22)
    plt.plot(gene_length, n1, label ="D1") 
    plt.plot(gene_length, n2, label = "D2")
    plt.plot(gene_length, n3, label = "D3")
    plt.plot(gene_length, n4, label = "D4") 
    plt.plot(gene_length, n5, label ="D5")
    plt.plot(gene_length, n6, label ="D6")
    plt.plot(gene_length, n7, label ="D7") 
    plt.plot(gene_length, n8, label ="D8")
    plt.legend(prop = { "size": 22 }, loc ="upper right")
    plt.savefig("gene_length.png", format="png")
    plt.clf()

if __name__ == "__main__":
    main()


