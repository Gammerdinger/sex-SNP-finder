Sex_SNP_finder
==============
Background
--------------
Sex_SNP_finder.pl was designed to find nucleotides that were fixed at a site in one population and polymorphic at the same site in a different population. It is also designed to look for regions of differentiation using a non-overlapping window approach. Initially, the program was used to compare between males and females, in order to find mutations potentially responsible for sex-determination. However, it is applicable to any population where a nucleotide is fixed in one pool and in a frequency of interest in another pool. Furthermore, functionality has been added to calculate Fst at each nucleotide that is polymorphic, not just those fitting the Sex_SNP_finder.pl criteria. This method was a derivation from Kofler et al. 2011's Popoolation2 calculation and treats deletions as a "5th nucleotide" when calculating Fst, a nuance not in Popoolation2. This calculation was not included in the BMC Genomics version of Sex_SNP_finder_now.pl

### Citation

If you are using the latest verion of Sex_SNP_finder (Sex_SNP_finder_GA.pl) please cite:

Gammerdinger WJ, Conte MA, Baroiller JF, D'Cotta H, Kocher TD: Comparative analysis of a sex chromosome from the blackchin tilapia, Sarotherodon melanotheron. BMC Genomics 2016 17(1):808.

Otherwise, cite Sex_SNP_finder as:

Gammerdinger WJ, Conte MA, Acquah EA, Roberts RB, Kocher TD: Structure and decay of a proto-Y region in Tilapia, Oreochromis niloticus. BMC Genomics 2014 15:975.

￼￼￼doi:10.1186/1471-2164-15-975

### Getting started with the Popoolation2 pipeline

This software is built off of the synchronized (.sync) file created in Popoolation2 (https://code.google.com/p/popoolation2/). Any issues getting to the .sync file should be able to be resolved there. An excellent tutorial on processing from .fastq files all the way through to the .sync file creation and further subsequent analyses (https://code.google.com/p/popoolation2/wiki/Tutorial) is provided there.

### Sex_SNP_finder.pl Features and Limitations

Features

* Able to find positions that are fixed or nearly fixed in one pool and within an assigned frequency in a different pool.
* Provides a non-overlapping window approach which can help identify regions that are the most dense of the SNPs that meet the parameters
* Provides a modified calculation of Fst from the one provided in Popoolation2. Popoolation 2 does not include deletions in its calculation of Fst like Sex_SNP_finder_Fst_now.pl does.

Limitations

* Popoolation2 provides a sliding window function for calculating Fst. Sex_SNP_finder_Fst_now.pl currently only calculates Fst for polymorphic bases and does not use windowing in its Fst calculation. 

### Command line example

```
perl Sex_SNP_finder_now.pl --input_file=test.sync --output_file=sex_SNP_finder_SNPs.igv --fixed_population=pool2 --fixed_threshold=0.9 --minimum_polymorphic_frequency=0.3 --maximum_polymorphic_frequency=0.7 --read_depth=10 --window_size=2 --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --description=description
```

### Command line input options

* input_file -> input_file.sync -> This is your sync file from the Popoolation2 pipeline. The format should look like this (You do not need to have the reference genome loaded into the 3rd column. Sex_SNP_Finder does not utilize reference genome information):
    - Column 1 -> Scaffold/Contig
    - Column 2 -> Position
    - Column 3 -> Reference genome base
    - Column 4 -> Pool1 A:T:C:G:N:Del
    - Column 5 -> Pool2 A:T:C:G:N:Del
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

* output_file -> output_file.igv -> This will be your IGV output file from Sex_SNP_finder.pl 
    - Column 1 -> Scaffold/Contig
    - Column 2 -> Start Position
    - Column 3 -> End Position
    - Column 4 -> Feature
    - Column 5 -> Frequency of SNP in non-fixed pool
        * We are looking for sites fixed in pool2 that are in an intermediate frequency (0.3-0.7) in pool1
        * Look at position 7 on scaffold_0
        * This position is fixed in pool2 for Thymine, but in pool1 1/3 of the positions have a Adenine and 2/3 have a Thymine
        * Therefore, the resulting IGV output looks like:

    
``` 
    Chromosome	Start	End	Feature	description_Sex_SNP_finder
    scaffold_0	7	8	snp	0.333333333333333
```

* non_overlapping_window_output_file -> This is an IGV-readable file where your non-overlapping window output will go. Following the parameters in the previous example along with paramters discussed below, the output should look like this:
   - Column 1 -> Scaffold
   - Column 2 -> Window start position
   - Column 3 -> Window end position
   - Column 3 -> Feature
   - Column 4 -> Number of Sex_SNP_finder.pl SNPs found in this window
      * Note: The first window is larger than the next three windows. This is because the first two rows of the .sync file fail to meet the read depth criteria.
```
   Chromosome   Start   End Feature description_Sex_SNP_finder
   scaffold_0	1	4	snp_count   0
   scaffold_0	5	6	snp_count   0
   scaffold_0	7	8	snp_count   1
   scaffold_0	9	10	snp_count   0
```

* fixed_population -> pool1 or pool2 -> This tells Sex_SNP_finder if you want it to consider pool 1 or pool2 to be the fixed pool.

* fixed_threshold -> Due to seqeuncing errors, you might want to consider not only sites that are fixed (those without a sequencing error) in your fixed population, but also sites that that are almost fixed (which may be the result of sequencing errors). We typically use a value of 0.9 here, but this is a function of coverage.

* minimum_polymorphic_frequency -> This is the lower bound on the SNP frequency in the polymorphic pool. In the previous example for the output file we used 0.3.

* maximum_polymoprhic_frequency -> This is the upper bound on the SNP frequency in the polymorphic pool. In the previous example for the output file we used 0.7.

* read_depth -> This is the minimum read depth that both pools must have before Sex_SNP_finder.pl looks for SNPs. We typically use a value of 10 here.

* window_size -> The size of your non-overlapping window. Sites with coverage that does not meet the read_depth threshold in both pools will not be considered. In the above example we used a window size of 2, but typically we use window sizes of 10,000. 

* description -> This is the description to be used as the header in the IGV file

### Sex_SNP_finder_Fst_now.pl

* The command line options are still same from Sex_SNP_finder_now.pl except now there is an fst_output_file option

### Command line example

```
perl Sex_SNP_finder_now.pl --input_file=test.sync --output_file=sex_SNP_finder_SNPs.igv --fst_output_file=fst.igv--fixed_population=pool2 --fixed_threshold=0.9 --minimum_polymorphic_frequency=0.3 --maximum_polymorphic_frequency=0.7 --read_depth=10 --window_size=2 --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --description=description
```

* New options

    - fst_output_file -> This is the IGV-readable Fst output file
        * Column 1 -> Scaffold
        * Column 2 -> Start Position
        * Column 3 -> End Position
        * Column 4 -> Feature
        * Column 5 -> Fst
```
Chromosome	Start	End	Feature	description_Fst
scaffold_0	3	4	snp	0.0595708639186899
scaffold_0	4	5	snp	0.160353535353535
scaffold_0	5	6	snp	0.173437908496732
scaffold_0	6	7	snp	0.138106995884774
scaffold_0	7	8	snp	0.212919254658385
scaffold_0	9	10	snp	0.0223166843783207
scaffold_0	10	11	snp	0.0554106783681708
```

Below is how Fst is calculated using Sex_SNP_finder_Fst_now.pl at each position:

![alt tag](https://github.com/Gammerdinger/sex-SNP-finder/blob/master/Popoolation2%20Fst%20equation.png)

### Features in the works

* Toggle options so that you can choose which computations you would like done.
* If you would like a feature added in, I am more than happy to collaborate and added features that you are interested in having.
