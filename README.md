sex_SNP_finder
==============
Background
--------------
Sex_SNP_finder.pl was designed to find nucleotides that were fixed at a site in one population and polymorphic at the same site in a different population. It is also designed to look for regions of differentiation using a sliding window approach. Initially, the program was used to compare between males and females, in order to find mutations potentially responsible for sex-determination. However, it is applicable to any population where one nucleotide is fixed in one pool and in a frequency of interest in another pool. Furthermore, functionality has been added to calculate Fst at each nucleotide that is polymorphic, not just those fitting the Sex_SNP_finder.pl criteria. This method was a derivation from Kofler et al. 2011's Popoolation2 calculation and treats deletions as a "5th nucleotide" when calculating Fst, a nuance not in Popoolation2. This calculation was not included in the BMC Genomics version of Sex_SNP_finder.pl

### Getting started with the Popoolation2 pipeline

This software is built off of the synchronized (.sync) file created in popoolation2 (https://code.google.com/p/popoolation2/). Any issues going from the the should be able to be resolved there. An excellent tutorial on processing from .fastq files all the way through the .sync file creation and further subsequent analyses (https://code.google.com/p/popoolation2/wiki/Tutorial) is provided there.

### Sex_SNP_finder.pl Features and Limitations

Features

-Able to find positions that are fixed or nearly fixed in one pool and within an assigned frequency in a different pool. -Provides a non-overlapping window approach which can help identify regions that are the most dense of the SNPs that meet the parameters -Provides a modified calculation of Fst from the one provided in Popoolation2. Popoolation 2 does not include deletions in its calculation of Fst.

Limitations

-Popoolation2 provides a sliding window function for calculating Fst. Sex_SNP_finder.pl currently only calculates Fst for polymorphic bases and does not use windowing in its Fst calculation. -Sex_SNP_finder.pl assumes a biallelic model. In other words, if your homogametic pool was fixed for adenine and your homogametic pool, with frequency thresholds of 0.3 to 0.7, had 0.5 adenine, 0.25 thymine and 0.25 guanine, this SNP would get called. Currently, we are addressing this situation, but in our data sets, this includes only a very small subset of sites called.

### Command line input options

* input_file -> input_file.sync -> This is your sync file from the Popoolation2 pipeline. The format should look like this (You do not need to have the reference genome loaded into the 3rd column. Sex_SNP_Finder.pl does utilize reference genome information):
    *Column 1 -> Scaffold/Contig
    *Column 2 -> Position
    *Column 3 -> Reference genome base
    *Column 4 -> Pool1 A:T:C:G:Del:N
    *Column 5 -> Pool2 A:T:C:G:Del:N
```
    scaffold_0	1	N	0:2:0:0:0:0	0:4:0:0:0:0
    scaffold_0	2	N	7:0:0:0:0:0	5:0:0:0:0:0
    scaffold_0	3	N	0:12:12:0:0:0	0:13:5:0:0:0
    scaffold_0	4	N	0:0:16:0:0:0	0:0:14:5:0:0
    scaffold_0	5	N	19:7:0:0:0:0	17:0:0:0:0:0
    scaffold_0	6	N	0:0:6:22:0:0	0:0:0:18:0:0
    scaffold_0	7	N	12:24:0:0:0:0	0:23:0:0:0:0
    scaffold_0	8	N	30:0:0:0:0:0	25:0:0:0:0:0
    scaffold_0	9	N	0:0:1:34:0:0	0:2:0:26:0:0
    scaffold_0	10	N	38:2:0:0:0:0	27:7:0:0:0:0
```

* --output_file -> output_file.igv -> This will be your IGV output file from Sex_SNP_finder.pl 
    *Column 1 -> Scaffold
    *Column 2 -> Start Position
    *Column 3 -> End Position
    *Column 4 -> Feature
    *Column 5 -> Frequency of SNP in non-fixed pool
        *We are looking for sites fixed in pool2 that are in an intermediate frequency (0.3-0.7) in pool1
        *Look at position 7 on scaffold_0
        *This position is fixed in pool2 for Thymine, but in pool1 1/3 of the positions have a Adenine and 2/3 have a Thymine

    
``` 
    scaffold_0	7	8	snp	0.333333333333333
```

* 
*
