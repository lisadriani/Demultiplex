#!/usr/bin/env python

import gzip
import bioinfo
import argparse
import numpy as np
from itertools import permutations


#make sure you're in py_310 module

def get_args():
    parser = argparse.ArgumentParser(description="Setting global variables")
    parser.add_argument("-R1",  "--readone", help="specify Read One", required = True)
    parser.add_argument("-R2", "--readtwo", help = "specify Read Two", required = True)
    parser.add_argument("-I1", "--indexone", help = "specify Index One", required = True)
    parser.add_argument("-I2", "--indextwo", help = "specify Index Two", required = True)
    parser.add_argument("-I", "--indexeslist", help = "Give a text file with the list of known indexes", required = True)
    return parser.parse_args()
	
args = get_args()
 
#Get all the known indexes and put them in a set called "known"
known = set()
revknown = set()
with open(args.indexeslist, "r") as index:
    for item in index: 
        item = item.strip()
        known.add(item)
for index in known: 
    revknown.add(bioinfo.revcomp(index))
#print(known)
#print(revknown)


#create a dictionary with the indexes as keys and the open file thing as the value 
fileR1 = {} #dictionary
fileR2 = {}
for item in known: 
    fileR1[item] = gzip.open(item+"_R1.fastq.gz", "wt")
    fileR2[item] = gzip.open(item+"_R2.fastq.gz", "wt")

#open all your input files bestie!
R1 = gzip.open(args.readone,"rt") #"t" needed so that it doesn't read the bits, but the actual letters/numbers
R2 = gzip.open(args.readtwo,"rt")
I1 = gzip.open(args.indexone,"rt")
I2 = gzip.open(args.indextwo,"rt")

#create and open the files for unknown and unmatched
unknownR1 = gzip.open("unknown_R1.fastq.gz", "wt") 
unknownR2 = gzip.open("unknown_R2.fastq.gz", "wt")
unmatchedR1 = gzip.open("unmatched_R1.fastq.gz", "wt")
unmatchedR2 = gzip.open("unmatched_R2.fastq.gz", "wt")

#create a dictionary to keep track of index pairs: 
#key = index_index possibility 
#value = number of times that happened
index_pairs = dict()
perms = permutations(known,2)
for x,y in perms:
    key = x+"_"+y
    index_pairs[key] = 0
#print(index_pairs)



#start a count for everything!
unmatch_count = 0
unknown_count = 0
low_qual = 0 
matched = dict()
for item in known:
    matched[item] = 0


i = 0
while True: 
    r1line1 = R1.readline().strip()
    r2line1 = R2.readline().strip()
    i1line1 = I1.readline().strip()
    i2line1 = I2.readline().strip()
    if r1line1 == "":
        break
    else:
        R1_array = np.array([r1line1,R1.readline().strip(),R1.readline().strip(),R1.readline().strip()])
        R2_array = np.array([r2line1,R2.readline().strip(),R2.readline().strip(),R2.readline().strip()])
        I1_array = np.array([i1line1,I1.readline().strip(),I1.readline().strip(),I1.readline().strip()])
        I2_array = np.array([i2line1,I2.readline().strip(),I2.readline().strip(),I2.readline().strip()])
        i +=1
        if i == 2:
            print("2!")
        if i == 10000000:
            print("Progress! 10,000,000 records done!")
            i = 0
        #now start checking the problems!


        # if "N" in I1_array[1]: 
        #     ##Write to unknown
        #     index1 = I1_array[1]
        #     revcomp = bioinfo.revcomp(I2_array[1])
        #     R1_array[0] = R1_array[0]+"_"+index1+"_"+revcomp+"\n"
        #     R2_array[0] = R2_array[0]+"_"+index1+"_"+revcomp+"\n"
        #     bioinfo.write_file(unknownR1,R1_array)
        #     bioinfo.write_file(unknownR2,R2_array)
        #     unknown_count += 1
        # if "N" not in I1_array[1]:
        #     if "N" in I2_array[1]:
        #         index1 = I1_array[1]
        #         revcomp = bioinfo.revcomp(I2_array[1])
        #         R1_array[0] = R1_array[0]+"_"+index1+"_"+revcomp+"\n"
        #         R2_array[0] = R2_array[0]+"_"+index1+"_"+revcomp+"\n"
        #         bioinfo.write_file(unknownR1,R1_array)
        #         bioinfo.write_file(unknownR2,R2_array)
        #         unknown_count += 1
            # else:
            #if "N" not in I2_array[1]:
        if I1_array[1] not in known:
            #index1 = I1_array[1]
            #revcomp = bioinfo.revcomp(I2_array[1])
            R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
            R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
            bioinfo.write_file(unknownR1,R1_array)
            bioinfo.write_file(unknownR2,R2_array)
            unknown_count += 1
            continue
            #write to unknown
        else:
        #if I1_array[1] in known:
            if I2_array[1] not in revknown:
                # index1 = I1_array[1]
                # revcomp = bioinfo.revcomp(I2_array[1])
                R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                bioinfo.write_file(unknownR1,R1_array)
                bioinfo.write_file(unknownR2,R2_array)
                unknown_count += 1
                continue
                #write to unknown
            else:
            #if bioinfo.revcomp(I2_array[1]) in known:
                if I1_array[1] != bioinfo.revcomp(I2_array[1]):
                    for index, value in enumerate(I1_array[3]):
                        Q1 = bioinfo.convert_phred(value)
                        Q2 = bioinfo.convert_phred(I2_array[3][index])
                        bool_break = False
                        if Q1 < 30 or Q2 <30:
                            R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                            R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                            bioinfo.write_file(unknownR1,R1_array)
                            bioinfo.write_file(unknownR2,R2_array)
                            bool_break = True
                            unknown_count +=1
                            low_qual += 1
                            break
                    if bool_break == True:
                        continue
                    else: 
                        index1 = I1_array[1]
                        revcomp = bioinfo.revcomp(I2_array[1])
                        R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                        R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                        bioinfo.write_file(unmatchedR1,R1_array)
                        bioinfo.write_file(unmatchedR2,R2_array)
                        index_pairs[index1+"_"+revcomp] +=1
                        unmatch_count += 1
                        continue
                else:
                #if I1_array[1] == bioinfo.revcomp(I2_array[1]):
                    for index, value in enumerate(I1_array[3]):
                        Q1 = bioinfo.convert_phred(value)
                        Q2 = bioinfo.convert_phred(I2_array[3][index])
                        bool_break = False
                        if Q1 < 30 or Q2 <30:
                            R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                            R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                            bioinfo.write_file(unknownR1,R1_array)
                            bioinfo.write_file(unknownR2,R2_array)
                            bool_break = True
                            unknown_count +=1
                            low_qual += 1
                            break
                    if bool_break == True:
                        continue
                    else:
                        for index, value in enumerate(I1_array[3]):
                            Q1 = bioinfo.convert_phred(value)
                            Q2 = bioinfo.convert_phred(I2_array[3][index])
                            bool_break = False
                            if Q1 < 30 or Q2 <30:
                                R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                                R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+bioinfo.revcomp(I2_array[1])
                                bioinfo.write_file(unknownR1,R1_array)
                                bioinfo.write_file(unknownR2,R2_array)
                                bool_break = True
                                unknown_count +=1
                                low_qual += 1
                                break
                        if bool_break == True:
                            continue
                        else: 
                            index1 = I1_array[1]
                            revcomp = bioinfo.revcomp(I2_array[1])
                            R1_array[0] = R1_array[0]+"_"+I1_array[1]+"_"+revcomp
                            R2_array[0] = R2_array[0]+"_"+I1_array[1]+"_"+revcomp
                            #print(fileR1[index1])
                            #print(fileR2[revcomp])
                            bioinfo.write_file(fileR1[index1],R1_array)
                            bioinfo.write_file(fileR2[revcomp],R2_array)
                            matched[index1] +=1
                            continue



R1.close()
R2.close()
I1.close()
I2.close()
unknownR1.close()
unknownR2.close()
unmatchedR1.close()
unmatchedR2.close()

for item in known: 
    fileR1[item].close()
    fileR2[item].close()

print("Number of Unknown records found:" , unknown_count)
print("Number of Unknown recods that are low quality:",low_qual)
print("Number of Unmatched Records found:" , unmatch_count)
print("Number of Matched Index Pairs by Index")
for key, value in matched.items():
    print(key,value)
print("Number of Hopped Indexes by Index Pair:")
for key, value in index_pairs.items():
    print(key,value) 







#print(fileR1)
#print(fileR2)
#dictionary with key as index, value as open()


    # print(itemr1)
    # print(itemr2)
    #create and open files with known file names for R1 and R2


# for item in known:
#     item+"_R1.fastq".write(R1_array[0]+"_"+index1+"_"+revcomp+"\n"+R1_array[1]+"\n"+R1_array[2]+"\n"+R1_array[3]+"\n")

#write to unknown

# R1_array = np.array([header,seq,+,Q])
# R2_array = np.array([header,seq,+,Q])
# I1_array = np.array([header,seq,+,Q])
# I2_array = np.array([header,seq,+,Q])

#don't open files in your loops! expensive and slow
# fh = open(args.file, "w")
# for i in range(10)
#     fh.write(str(i) + leslie is awesome)

        # print(R1_array)
        # print(R2_array)
        # print(I1_array)
        # print(I2_array)
        # print(I2_array[1])
        # print(I1_array[1])
        # print(I2_array[2])
        # print(I2_array[3])

        # r1line2 = R1.readline()
        # r1line2 = r1line2.strip()
        # r1line3 = R1.readline()
        # r1line3 = r1line3.strip()
        # r1line4 = R1.readline()
        # r1line4 = r1line4.strip()
                # r2line2 = R2.readline()
        # r2line2 = r2line2.strip()
        # r2line3 = R2.readline()
        # r2line3 = r2line3.strip()
        # r2line4 = R2.readline()
        # r2line4 = r2line4.strip()

        # R2_array = np.array([r2line1,r2line2,r2line3,r2line4])
        # i1line2 = I1.readline()
        # i1line2 = i1line2.strip()
        # i1line3 = I1.readline()
        # i1line3 = i1line3.strip()
        # i1line4 = I1.readline()
        # i1line4 = i1line4.strip()
        # I1_array = np.array([i1line1,i1line2,i1line3,i1line4])
        # i2line2 = I2.readline()
        # i2line2 = i2line2.strip()
        # i2line3 = I2.readline()
        # i2line3 = i2line3.strip()
        # i2line4 = I2.readline()
        # i2line4 = i2line4.strip()
        # I2_array = np.array([i2line1,i2line2,i2line3,i2line4])

# fileR1 = index1+"_R1.fastq"
# fileR2 = revcomp+"_R2.fastq"


# unknownR1.write(R1_array[0]+R1_array[1]+"\n"+R1_array[2]+"\n"+R1_array[3]+"\n")
# unknownR2.write(R2_array[0]+R2_array[1]+"\n"+R2_array[2]+"\n"+R2_array[3]+"\n")

#unmatchedR1.write(R1_array[0]+R1_array[1]+"\n"+R1_array[2]+"\n"+R1_array[3]+"\n")
# unmatchedR2.write(R2_array[0]+R2_array[1]+"\n"+R2_array[2]+"\n"+R2_array[3]+"\n")

#fileR1[index1].write(R1_array[0]+R1_array[1]+"\n"+R1_array[2]+"\n"+R1_array[3]+"\n")
# fileR2[revcomp].write(R1_array[0]+R1_array[1]+"\n"+R1_array[2]+"\n"+R1_array[3]+"\n")

