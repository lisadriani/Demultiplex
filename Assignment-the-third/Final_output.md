After demultiplexing reads from 2017 sequencing data: 

226,715,602 indexes were within the threshold >30 Q score, no N's, within the list of known indexes, and showed paired end matching. 

Index swapping happened 330,975 times. 
Of the 136,200,158 "unknown" records, 105,416,196 (77.4%) of them were due to low quality index reads. 


## Number of Matched Index Pairs by Index
### Percentage of each index pair 

```
Index -- Count -- % of total Matched Index Pairs
CGATCGAT 4237854 1.87%
AGGATAGC 5861709 2.59%
AACAGCGA 6368144 2.81%
TCGAGAGT 7448072 3.29%
ATCATGCG 6927867 3.06%
GCTACTCT 4301318 1.90%
AGAGTCCA 7602663 3.35%
ACGATCAG 5933528 2.62%
TCTTCGAC 30089661 13.27%
TAGCCATG 7148153 3.15%
TACCGGAT 49686878 21.92%
ATCGTGGT 4730009 2.09%
GTCCTAAG 6200133 2.73%
TCGACAAG 2644260 1.17%
CTAGCTCA 13034311 5.75%
CTCTGGAT 24515042 10.81%
TGTTCCGT 11450554 5.05%
CACTTCAC 2577666 1.14%
TCGGATTC 2874320 1.27%
GATCAAGG 4628196 2.04%
CGGTAATC 2393021 1.06%
TATGGCAC 7651472 3.37%
GTAGCGTA 5774439 2.55%
GATCTTGC 2636332 1.16%

```

## Distribution of Index reads

!["Assignment-the-third/pieplot.png"](https://github.com/lisadriani/Demultiplex/blob/master/Assignment-the-third/pieplot.png)


## Count of matched index reads

![bargraph](https://github.com/lisadriani/Demultiplex/blob/master/Assignment-the-third/pairedd.png)


## Number of Hopped Indexes by Index Pair:
### Percentages are calculated with a total of all hopped indexes
**Only lists those above 1%**

```
Index Pair  --  Count -- Percentage across all hopped indexes

TATGGCAC_TGTTCCGT 58741 17.75%
TGTTCCGT_TATGGCAC 58569 17.70%
GATCAAGG_TCTTCGAC 8799 2.66%
CTAGCTCA_TCGACAAG 8656 2.62%
TCGACAAG_ATCATGCG 5838 1.76%
TACCGGAT_CTCTGGAT 5609 1.69%
TACCGGAT_TCTTCGAC 5563 1.68%
CTCTGGAT_TACCGGAT 5463 1.65%
GTCCTAAG_TATGGCAC 4557 1.38%
TCTTCGAC_TACCGGAT 4133 1.25%
CGGTAATC_TACCGGAT 3614 1.09%

```
