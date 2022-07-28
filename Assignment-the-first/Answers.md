# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 | 101 | Phred 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 | 8 | Phred 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 8 | Phred 33 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 | 101 | Phred 33 |

2. Per-base NT distribution
    1. ![read2](https://user-images.githubusercontent.com/70485602/181601893-0cd6748c-49b4-44d8-be24-a2a91d308ffb.png)
![read1](https://user-images.githubusercontent.com/70485602/181601914-85f10e7a-7b4b-4e64-967d-a95c8caff659.png)
![index2](https://user-images.githubusercontent.com/70485602/181601934-58cf008c-19a1-4f93-82ae-ac6b850577a5.png)
![index1](https://user-images.githubusercontent.com/70485602/181601960-85c4790c-cda8-4555-bd18-561a972e85af.png)

    2. An individual Qscore of 30 for the index reads are important for sample identification because it is important to truly get the sample that is associated with your experiment. If there is more than a chance than a 1 in 1000 chance that the nucelotide is incorrect in a sequence, that could lead to major downstream analysis issues especially if you are sequencing a new genome. For the sequence itself, I think the Qscore cutoff can be a little bit lower, but still ideally Q>30 because most genome assembly programs have some wiggle room to allow for sequences that ~almost~ match the known sequences. Also, if you have a read pair, then you can compare the reverse complement of that strand to double check the sequence! 
    3. 3976613 in Index 1
       3328051 in Index 2
    
## Part 2
1. Define the problem
When multiplexing sequences, the sequences have barcodes on either end of the DNA sequences that will allow the sequences to be sorted after to where it belongs. However, during the reaction, these barcodes could incorrectly bind and a barcode that doesn't match the sequence could bind to the sequence. These samples specifically have the same barcode on either side, so if the index/barcodes don't match on either end of the sequence, we know they were a victim of index hopping. Read 1 and Read 2 are the reverse complements of each other and Index 1 and Index 2 are also reverse compelements. 
2. Describe output
The output will have 24 Read 1 matched files, 24 Read 2 matched files, 1 Read 1 unmatched (known) file, 1 Read 2 unmatched (known) file, 1 Read 1 unknown file, 1 read unknown file. The 24 matched files correspond to known indexes. 
"matched" means that index 1 and index 2 are matched to each other (and known in the our list of 24 indexes), so we know there was no index hopping
"unmatched" means that index 1 and index 2 are in our known list of 24 indexes, but the index 1 is not the same as index 2. 
"unknown" means that the qscore is lower than 30 and/or one/both of the indexes are not in the known list of 24 (which also means there could be an "N")
Each file will hold the entire record of the reads, with an edited header that has "header"+_index1_revcomp(index2). I want the reverse compelement of read 2 in the header so that when I'm using the matched files downstream, I know that the indexes are matched and which indexes there are. This will also make it easier to know exactly which file I'm looking at just by looking at the first few records in each of the files. 
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement


R1_array = [header,seq,+,Q]

set = {known indexes}
unknown_count = 0
unmatched_count = 0
match_count = 0
Open the four fastq files and read each simultaenously
    for each nth record of each file--> each file has an array size of 4 for each line #take each /record/ at a time
        check that phred score for each position in index 1 and index 2 are <30 
            if not: write to r1 and r2 unknown files with edited headers
                edited headers = index1_revcomp(index2)
                unknown_count+=1
        check index 1 in set?
            if not: send to unknown
                edited headers = index1_revcomp(index2)
                unknown_count+=1
            if yes: is revcomp(index 2) in set?  #set rev comp to a variable to write it
                if not: send to unknown
                    edited headers = index1_revcomp(index2)
                    unknown_count+=1
                if yes: check index1 to revcomp(index2) matching
                    if not: send to unmatch
                        edited headers = index1_revcomp(index2)
                        unmatch_count += 1
                    if yes: write to individual matched file (labeled with index1_revcomp(index2)?)
                        edited headers = index1_revcomp(index2) 
                        match_count +=1 






Functions: 

def revcomp(seq: str): 
    ```takes a sequence string and returns the reverse complement``` 
    return revcomp(seq)

Test example: revcomp(AAACCCTTTGGG) = CCCAAAGGGTTT

def convert_phred(letter: str) -> int:
    ```Converts a single character into a phred score```
    return (ord(letter)) - 33

Test example: convert_phred(E) = 36


