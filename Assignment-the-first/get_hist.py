#!/usr/bin/env python

# /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz


import numpy as np
import bioinfo
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Setting global variables")
    parser.add_argument("-l", "--readlength", help = "specify read length of sequence line", required = True, type=int)
    parser.add_argument("-f",  "--filename", help="specify input filename", required = True)
    parser.add_argument("-o", "--output", help = "specify output file name for histogram", required = True)
    return parser.parse_args()
	
args = get_args()

# num_lines = 0
# with gzip.open(args.filename, "r") as f:
#     for line in f:
#         num_lines+=1
# all_qscores = np.zeros((101,num_lines//4), dtype=int)

# all_qscores = [] * args.readlength 
# my_list: list = []
# my_list = bioinfo.init_list(my_list, args.readlength)
#mean = np.zeros(101,dtype=np.float)

def populate_list(file: str) -> tuple[list, int]:
    """This function will convert each quality score into it's corresponding number and add it to an ongoing quality 
    score for each sequence position, and return the list of summed quality scores and the number of lines in the file"""
    with gzip.open(args.filename,"r") as f:
        my_list = [0] * args.readlength
        #my_list = bioinfo.init_list(my_list, args.readlength)
        i = 0 
        for line in f: 
            i += 1      
            if i %4 == 0:
                line = line.strip()
                line = line.decode("ascii")
                for count,letter in enumerate(line):
                    converted = bioinfo.convert_phred(letter)
                    my_list[count] += converted
    return my_list, i

my_list, num_lines = populate_list(args.filename)

# xlist = populate_list(args.filename)
# zlist = xlist[0]
for index,line in enumerate(my_list):
        my_list[index] = line / (num_lines//4)
        #print(f"{num_lines}\t{zlist[num_lines]}")

# with gzip.open(args.filename, "r") as f:
#     i = 0
#     for line in f:
#         i+=1
#         if i %4 == 0:
#             line = line.strip()
#             line = line.decode("ascii")
#             counter_array = i//4
#             for count,letter in enumerate(line):
#                 converted = bioinfo.convert_phred(letter) #converting each score into a qscores
#                 all_qscores[count, counter_array -1 ] = converted #filling the qscores array

#print(zlist[num_lines])
# with open(file,"r") as f:
#     i = 0 
#     for line in f: 
#         i += 1      
#         if i %4 == 0:
#             line = line.strip()
#             counter_array = i//4
#             for count,letter in enumerate(line):
#                 converted = bioinfo.convert_phred(letter) #converting each score into a qscores
#                 all_qscores[count, counter_array-1] = converted #filling the qscores array

# for index, element in enumerate(all_qscores):
#     mean[index] = np.mean(element)

#print(mean)
import matplotlib.pyplot as plt

x = range(int(args.readlength))
y = my_list

fig, ax = plt.subplots()
ax.bar(x,y)
plt.xlabel("Position")
plt.ylabel("Average Q-Score")
plt.title("Average Q-Score per position")

plt.savefig(str(args.output)+".png") #save the file

