sex_SNP_finder
==============
Background
--------------
Sex_SNP_finder.pl was designed to find nucleotides that were fixed at a site in one population and polymorphic at the same site in a different population. It is also designed to look for regions of differentiation using a sliding window approach. Initially, the program was used to compare between males and females, in order to find mutations potentially responsible for sex-determination. However, it is applicable to any population where one nucleotide is fixed in one population and polymorphic in the other. Furthermore, functionality has been added to calculate Fst at each nucleotide that is polymorphic, not just those fitting the Sex_SNP_finder.pl criteria. This method was a derivation from Kofler et al. 2011's Popoolation2 calculation and treats deletions as a "5th nucleotide" when calculating Fst, a nuance not in Popoolation2.
