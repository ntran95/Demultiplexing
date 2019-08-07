#!/usr/bin/env python

#import numpy
import argparse
import gzip

def get_arguements():
    parser = argparse.ArgumentParser(description="reading in different files")
    parser.add_argument("-f", "--file_name", help="this argument specifies the filename", type =str, required=True)
    parser.add_argument("-l", "--bp_length", help="this argument specifies the base pair length in each seq", type =int, required=True)
    parser.add_argument("-p", "--plot", help="this argument specifies the title of the plot", type =str, required=True)
    return parser.parse_args()

args=get_arguements()
f=args.file_name
l=args.bp_length
p=args.plot

def convert_phred(letter):
    """Converts a single character into a phred score"""
    return (ord(letter)) - 33

mean_scores=[]
for i_mean in range(l):
    mean_scores.append(i_mean)
print(len(mean_scores))

'''open the FASTQ file and loop through every record (recommend testing your code with a smaller subsample of the file). Convert the Phred quality score from a letter to its corresponding number and add it to an ongoing sum of the quality scores for each base pair. So, the quality score of the first nucleotide of every read will be summed together in position 0 of the array you create. Likewise, the quality scores of the 101th nucleotide will be stored in position 100 of the array'''
#with open(f, "r") as fh:
with gzip.open(f,"rt") as fh:
    LN = 0
    for line in fh:
        line = line.strip('\n')
        LN+=1
        if LN%4 == 0:
            index_location=0
                #print(line)
            for i in line:
                results = convert_phred(i)
                #print(results)
                #sums up quality scores between the base pair locations
                mean_scores[index_location] = mean_scores[index_location] + results
                index_location+=1


index_location = 0
print ("# of Base Pair" "\t" "Mean Quality Score")
for i in mean_scores:
    mean_scores[index_location] = mean_scores[index_location]/(LN/4)
    print(index_location,"\t", mean_scores[index_location])
    index_location+=1

import matplotlib.pyplot as plt
y= mean_scores
x = range(l)
plt.bar(x,y)
plt.ylabel('Mean Scores')
plt.xlabel('# of base pairs')
#plt.show()
plt.savefig(p)
plt.close()
#plt.close(fig)
print("finished plot")
