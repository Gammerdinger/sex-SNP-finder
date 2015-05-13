#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(min);
use Math::Complex;

=pod
 
=head1 NAME
 
Sex_SNP_finder_GA.pl
 
=head1 AUTHORS
 
 Will Gammerdinger - Program Designer and Matt Conte - Coding Advisor
 
=head1 DESCRIPTION
 
 Sex_SNP_finder_GA.pl was a derivation from Sex_SNP_finder.pl. It is designed to look for regions of differentiation using a non-overlapping window approach. Sex_finder.pl was designed to find nucleotides that were fixed at a site in one population and polymorphic at the same site in a different population. Initially, the program was used to compare between males and females, in order to find mutations responsible for sex-determination. However, it is applicable to any population where one nucleotide is fixed in one population and polymorphic in the other.
 
=head1 EXAMPLE
 
 perl Sex_SNP_finder_GA.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --minimum_read_depth=[value greater than 0] --maximum_read_depth=[value greater than the minimum_read_depth] --minimum_allele_count=[minimum number of allele counts at a site before an allele is considered real and not a sequencing error] --sex_SNP_finder_window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --fst_output_file=fst_output_file.igv --dxy_output_file=dxy_output_file.igv --da_output_file=da_output_file.igv --cp_output_file=cp_output_file.igv --neis_D_output_file=neis_D_output_file.igv--description=description of file to be used in IGV header [--help|-?]

 
=head1 LICENSE
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
=head1 VERSION
 
 version0.1.0 - Sex_SNP_finder_Fst_Dxy.pl was updated to Sex_SNP_finder_GA.pl with the GA standing for Genomics Analysis. Several key features were added. First, a minimum allele count feature was added in order to try to filter out sequencing errors. Setting this to 0 or 1 will therefore not filter out any low frequency SNPs. Additionally, a maximum read depth feature was added to exclude sites with high coverage and also allow for an upper bound on Nei's D to be calculated. Calculations for Da, Cp and Nei's D have been added. Formulas for these equations can be found on the GitHub formula section. Nei's D for fixed sites is infinity, which isn't terrible helpful when it is trying to be viewed in a genome browser such as IGV. So, a "maximum Nei's D value" is calculated based upon one population being fixed at the maximum read depth and the other population being almost alternatively fixed for the opposite allele but it has one count in common with the former population. For example if the maximum_read_depth is equal to 100 then the "maximum Nei's D value" would be calculated for population 1 being fixed for 100 Adenine counts and population 2 would be polymorphic for 99 Thymine counts and 1 Adenine count. In this example the "maximum Nei's D value" is equal to 4.59517.
 
 version0.0.8 - Dxy funcationality was added and now sites with scores of 0 for Fst and Dxy are printed. They can be easily filtered out in downstream analysis using code on the Github page.
 
 version0.0.7 - The biallelic assumption was corrected. Additionally, the non-overlapping window output was changed to an IGV readable output.
 
 version0.0.6 - Fst component was added at the request of Matt Conte. It was calculated in a manner similar to Kofler et al., 2011b. However, Kofler et al., 2011 did not include frequency of deletions in their Fst calculations. We have made this alteration. The name was changed to Sex_SNP_finder_Fst_now.pl to reflect this change.
 
 version0.0.5 - Names of Sex_finder software were changed from Sex_finder to Sex_SNP_finder. Added a non-overlapping window approach to Sex_SNP_finder.pl to make Sex_SNP_finder_now.pl. This non-overlapping window skips sites that fail to meet the read depth threshold and gaps.
 
 version0.0.4 - Returned a guide prompt when option fields were left empty.
 
 version0.0.3 - Changed the input and output from arguments to options.
 
 version0.0.2 - Allowed declaration of pooled data set
 
 version0.0.1 - Original Version
 
=cut


my ($input_file, $output_file, $fixed_population, $fixed_threshold, $minimum_polymorphic_frequency, $maximum_polymorphic_frequency, $minimum_read_depth, $maximum_read_depth, $minimum_allele_count,  $sex_SNP_finder_window_size, $non_overlapping_window_output_file, $fst_output_file, $dxy_output_file, $da_output_file, $cp_output_file, $neis_D_output_file, $description, $help) = ("empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty","empty", "empty", "empty", "empty", "empty");

GetOptions(
"input_file=s"                              => \$input_file,
"output_file=s"                             => \$output_file,
"fixed_population=s"                        => \$fixed_population,
"fixed_threshold=s"                         => \$fixed_threshold,
"minimum_polymorphic_frequency=s"           => \$minimum_polymorphic_frequency,
"maximum_polymorphic_frequency=s"           => \$maximum_polymorphic_frequency,
"minimum_read_depth=s"                      => \$minimum_read_depth,
"maximum_read_depth=s"                      => \$maximum_read_depth,
"minimum_allele_count=s"                    => \$minimum_allele_count,
"sex_SNP_finder_window_size=s"              => \$sex_SNP_finder_window_size,
"non_overlapping_window_output_file=s"      => \$non_overlapping_window_output_file,
"fst_output_file=s"                         => \$fst_output_file,
"dxy_output_file=s"                         => \$dxy_output_file,
"da_output_file=s"                          => \$da_output_file,
"cp_output_file=s"                          => \$cp_output_file,
"neis_D_output_file=s"                      => \$neis_D_output_file,
"description=s"                             => \$description,
"help|?"                                    => \$help
) or Usage ( "Invalid command-line option.");

Usage() if defined $help;

if ($input_file eq "empty" || $output_file eq "empty" || $fixed_population eq "empty" || $fixed_threshold eq "empty" || $minimum_polymorphic_frequency eq "empty" || $maximum_polymorphic_frequency eq "empty" || $minimum_read_depth eq "empty" || $maximum_read_depth eq "empty" || $minimum_allele_count eq "empty" ||   $sex_SNP_finder_window_size eq "empty" || $non_overlapping_window_output_file eq "empty" || $fst_output_file eq "empty" ||  $dxy_output_file eq "empty" || $da_output_file eq "empty" || $cp_output_file eq "empty" || $neis_D_output_file eq "empty" || $description eq "empty"){
    die "\nERROR: The format should be perl Sex_SNP_finder_GA.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --minimum_read_depth=[value greater than 0] --maximum_read_depth=[value greater than the minimum_read_depth] --minimum_allele_count=[minimum number of allele counts at a site before an allele is considered real and not a sequencing error] --sex_SNP_finder_window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --fst_output_file=fst_output_file.igv --dxy_output_file=dxy_output_file.igv --da_output_file=da_output_file.igv --cp_output_file=cp_output_file.igv --neis_D_output_file=neis_D_output_file.igv--description=description of file to be used in IGV header [--help|-?]\n\nOne or more of your option fields is empty.\n\nThe parameters are currently:\n\ninput_file: $input_file\noutput_file: $output_file\nfixed_population: $fixed_population\nfixed_threshold: $fixed_threshold\nminimum_polymorphic: $minimum_polymorphic_frequency\nmaximum_polymorphic_frequency: $maximum_polymorphic_frequency\nminimum_read_depth: $minimum_read_depth\nmaximum_read_depth: $maximum_read_depth\nminimum_allele_count: $minimum_allele_count\nsex_SNP_finder_window_size: $sex_SNP_finder_window_size\nnon_overlapping_window_output_file: $non_overlapping_window_output_file\nfst_output_file: $fst_output_file\ndxy_output_file: $dxy_output_file\nda_output_file: $da_output_file\ncp_output_file: $cp_output_file\nneis_D_output_file: $neis_D_output_file\ndescription: $description\n\nFor more information, use the command perldoc Sex_SNP_finder_Fst_now.pl\n\n"
}

open (my $INPUT, "<$input_file");

open (my $OUTPUT, ">$output_file");
open (my $SL_OUTPUT, ">$non_overlapping_window_output_file");
open (my $FST_OUTPUT, ">$fst_output_file");
open (my $DXY_OUTPUT, ">$dxy_output_file");
open (my $DA_OUTPUT, ">$da_output_file");
open (my $CP_OUTPUT, ">$cp_output_file");
open (my $NEIS_D_OUTPUT, ">$neis_D_output_file");


print "Your fixed population is ", $fixed_population, ".\n";
print "Your minimum read depth is ", $minimum_read_depth, " nucleotides.\n";
print "Your maximum read depth is ", $maximum_read_depth, " nucleotides.\n";
print "Your minimum allele count is ", $minimum_allele_count, " counts.\n";
print "Your minimum fixed threshold is ", $fixed_threshold, ".\n";
print "Your minumum polymorphic frequency is ", $minimum_polymorphic_frequency, ".\n";
print "Your maximum polymorphic frequency is ", $maximum_polymorphic_frequency, ".\n";
print "Your sex_SNP_finder non-overlapping window is ", $sex_SNP_finder_window_size, ".\n";


if ($minimum_polymorphic_frequency > $maximum_polymorphic_frequency){
    die "ERROR: The minimum polymorphic frequency is greater than the maximum polymorphic frequency.\nRemember that these allele frequencies of all sites that are not the same as the fixed allele.\n"
}

if ($minimum_read_depth > $maximum_read_depth){
    die "ERROR: The minimum read depth is greater than the maximum read depth."
}

my $first_linkage_group;
my $first_position;
my $last_linkage_group;
my $last_position;
my $boolean_1 = "false";
my $boolean_2 = "false";

while (my $line = <$INPUT>){
    $line =~ s/:/\t/g;
    my @array_of_line = split(/\t/, $line);
    
    my $Scaffold = $array_of_line[0];
    my $Start_position = $array_of_line[1];
    if ($boolean_1 eq "false"){
        $first_linkage_group = $array_of_line[0];
        $boolean_1 = "true"
    }
    if ($boolean_2 eq "false"){
        $first_position = $array_of_line[1];
        $boolean_2 = "true"
    }
    $last_linkage_group = $array_of_line[0];
    $last_position = $array_of_line[1]
}

print "Your input file begins at ", $first_linkage_group, " ", $first_position, " and ends at ", $last_linkage_group," ", $last_position, ".\n";

close $INPUT;

open ($INPUT, "<$input_file");

print $OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Sex_SNP_finder\n";

print $SL_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Sex_SNP_finder_count\n";

print $FST_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Fst\n";

print $DXY_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Dxy\n";

print $DA_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Da\n";

print $CP_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Cp\n";

print $NEIS_D_OUTPUT "Chromosome\tStart\tEnd\tFeature\t$description"."_Nei's_D\n";

my $i = 0;
my $window_counter=0;
my $snp_counter=0;
my $old_scaffold = $first_linkage_group;
my $start_position_boolean = "false";
my $window_first_position;
my $window_last_position;

my $Adenine_count_fixed_max  = $maximum_read_depth - 1;
my $Thymine_count_fixed_max  = 1;
my $Cytosine_count_fixed_max = 0;
my $Guanine_count_fixed_max  = 0;
my $Deletion_count_fixed_max = 0;
my $Adenine_count_polymorphic_max  = 0;
my $Thymine_count_polymorphic_max  = $maximum_read_depth;
my $Cytosine_count_polymorphic_max = 0;
my $Guanine_count_polymorphic_max  = 0;
my $Deletion_count_polymorphic_max = 0;
my $read_nucleotides_fixed_max = $maximum_read_depth;
my $read_nucleotides_polymorphic_max = $maximum_read_depth;
my $neis_D_max = 0;

my ($Dxy_value_ignore, $Cp_value_ignore, $maximum_neis_D) = Dxy($Adenine_count_fixed_max, $Thymine_count_fixed_max, $Cytosine_count_fixed_max, $Guanine_count_fixed_max, $Deletion_count_fixed_max, $read_nucleotides_fixed_max, $Adenine_count_polymorphic_max, $Thymine_count_polymorphic_max, $Cytosine_count_polymorphic_max, $Guanine_count_polymorphic_max, $Deletion_count_polymorphic_max, $read_nucleotides_polymorphic_max, $neis_D_max);


if ($fixed_population =~ m/pool1/){
    while (my $line = <$INPUT>){
        $line =~ s/:/\t/g;
        my @array_of_line = split(/\t/, $line);
        my $Scaffold = $array_of_line[0];
        my $Start_position = $array_of_line[1];
        my $End_position   = $array_of_line[1] + 1;
        my $Feature = "snp";
        my $SL_feature = "snp_count";
        my $Adenine_count_fixed  = $array_of_line[3];
        my $Thymine_count_fixed  = $array_of_line[4];
        my $Cytosine_count_fixed = $array_of_line[5];
        my $Guanine_count_fixed  = $array_of_line[6];
        my $Deletion_count_fixed = $array_of_line[8];
        my $Adenine_count_polymorphic  = $array_of_line[9];
        my $Thymine_count_polymorphic  = $array_of_line[10];
        my $Cytosine_count_polymorphic = $array_of_line[11];
        my $Guanine_count_polymorphic  = $array_of_line[12];
        my $Deletion_count_polymorphic = $array_of_line[14];
        my $Fst_boolean = 0;
        my $Fst_value = 0;
        my $Dxy_boolean = 0;
        my $Dxy_value = 0;
        my $Pi_fixed = 0;
        my $Pi_polymorphic = 0;
        my $Da_value = 0;
        my $Cp_value = 0;
        my $Neis_D_value = 0;
        
        if ($Adenine_count_fixed <= $minimum_allele_count){
            $Adenine_count_fixed = 0
        }
        if ($Thymine_count_fixed <= $minimum_allele_count){
            $Thymine_count_fixed = 0
        }
        if ($Cytosine_count_fixed <= $minimum_allele_count){
            $Cytosine_count_fixed = 0
        }
        if ($Guanine_count_fixed <= $minimum_allele_count){
            $Guanine_count_fixed = 0
        }
        if ($Deletion_count_fixed <= $minimum_allele_count){
            $Deletion_count_fixed = 0
        }
        if ($Adenine_count_polymorphic <= $minimum_allele_count){
            $Adenine_count_polymorphic = 0
        }
        if ($Thymine_count_polymorphic <= $minimum_allele_count){
            $Thymine_count_polymorphic = 0
        }
        if ($Cytosine_count_polymorphic <= $minimum_allele_count){
            $Cytosine_count_polymorphic = 0
        }
        if ($Guanine_count_polymorphic <= $minimum_allele_count){
            $Guanine_count_polymorphic = 0
        }
        if ($Deletion_count_polymorphic <= $minimum_allele_count){
            $Deletion_count_polymorphic = 0
        }
        
        my $read_nucleotides_fixed = $Adenine_count_fixed + $Thymine_count_fixed + $Cytosine_count_fixed + $Guanine_count_fixed + $Deletion_count_fixed;
        my $read_nucleotides_polymorphic = $Adenine_count_polymorphic + $Thymine_count_polymorphic + $Cytosine_count_polymorphic + $Guanine_count_polymorphic + $Deletion_count_polymorphic;
        
        if ($read_nucleotides_fixed > $maximum_read_depth || $read_nucleotides_polymorphic > $maximum_read_depth){
            next;
        }
        
        if ($read_nucleotides_fixed < $minimum_read_depth || $read_nucleotides_polymorphic < $minimum_read_depth){
            next;
        }
        
        if ($start_position_boolean eq "false"){
            $window_first_position = $array_of_line[1];
            $start_position_boolean = "true"
        }
        
        if ($Scaffold ne "$old_scaffold"){
            if ($window_counter!= 0 && $window_counter >= 0.1 * $sex_SNP_finder_window_size){
                $snp_counter = $snp_counter / $window_counter * $sex_SNP_finder_window_size;
                print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n"
            }
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
        
        $old_scaffold = $Scaffold;
        $window_last_position = $array_of_line[1];
        
        
        if ($read_nucleotides_fixed >= $minimum_read_depth){
            if($read_nucleotides_polymorphic >= $minimum_read_depth){
                # window counter
                $window_counter = $window_counter + 1;
                
                # Adenine
                
                my $Adenine_proportion_fixed = $Adenine_count_fixed / $read_nucleotides_fixed;
                my $Non_adenine_proportion_polymorphic = 1 - $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_adenine_proportion_polymorphic > 0 && $Adenine_proportion_fixed <= 1){
                    if($Non_adenine_proportion_polymorphic <= 1 && $Adenine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Adenine_proportion_fixed >= $fixed_threshold){
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                    }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                    }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                    }
                }
                
                
                #Thymine
                
                my $Thymine_proportion_fixed = $Thymine_count_fixed / $read_nucleotides_fixed;
                my $Non_thymine_proportion_polymorphic = 1 - $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_thymine_proportion_polymorphic > 0 && $Thymine_proportion_fixed <= 1){
                    if($Non_thymine_proportion_polymorphic <= 1 && $Thymine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Thymine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                        }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                    }
                }
                
                # Cytosine
                
                my $Cytosine_proportion_fixed = $Cytosine_count_fixed / $read_nucleotides_fixed;
                my $Non_cytosine_proportion_polymorphic = 1 - $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_cytosine_proportion_polymorphic > 0 && $Cytosine_proportion_fixed <= 1){
                    if($Non_cytosine_proportion_polymorphic <= 1 && $Cytosine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Cytosine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                    }
                }
                
                # Guanine
                
                my $Guanine_proportion_fixed = $Guanine_count_fixed / $read_nucleotides_fixed;
                my $Non_guanine_proportion_polymorphic = 1 - $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_guanine_proportion_polymorphic > 0 && $Guanine_proportion_fixed <= 1 ){
                    if ($Non_guanine_proportion_polymorphic <= 1 && $Guanine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Guanine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                    }
                }
                
                # Deletion
                
                my $Deletion_proportion_fixed = $Deletion_count_fixed / $read_nucleotides_fixed;
                my $Non_deletion_proportion_polymorphic = 1 - $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_deletion_proportion_polymorphic > 0 && $Deletion_proportion_fixed <= 1){
                    if ($Non_deletion_proportion_polymorphic <= 1 && $Deletion_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                            $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Deletion_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                    }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                }
            }
        }
        print $FST_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Fst_value\n";
        print $DXY_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Dxy_value\n";
        print $DA_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Da_value\n";
        print $CP_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cp_value\n";
        print $NEIS_D_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Neis_D_value\n";
        $i = $i + 1;
        my $d = 1000000;
        if ($i % $d == 0){
            print "[$i] nucleotide positions analyzed.\n";
        }
        if ($window_counter == $sex_SNP_finder_window_size){
            print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n";
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
    }
}

elsif ($fixed_population =~ m/pool2/){
    while (my $line = <$INPUT>){
        $line =~ s/:/\t/g;
        my @array_of_line = split(/\t/, $line);
        my $Scaffold = $array_of_line[0];
        my $Start_position = $array_of_line[1];
        my $End_position   = $array_of_line[1] + 1;
        my $Feature = "snp";
        my $SL_feature = "snp_count";
        my $Adenine_count_fixed  = $array_of_line[9];
        my $Thymine_count_fixed  = $array_of_line[10];
        my $Cytosine_count_fixed = $array_of_line[11];
        my $Guanine_count_fixed  = $array_of_line[12];
        my $Deletion_count_fixed = $array_of_line[14];
        my $Adenine_count_polymorphic  = $array_of_line[3];
        my $Thymine_count_polymorphic  = $array_of_line[4];
        my $Cytosine_count_polymorphic = $array_of_line[5];
        my $Guanine_count_polymorphic  = $array_of_line[6];
        my $Deletion_count_polymorphic = $array_of_line[8];
        my $Fst_boolean = 0;
        my $Fst_value = 0;
        my $Dxy_boolean = 0;
        my $Dxy_value = 0;
        my $Pi_fixed = 0;
        my $Pi_polymorphic = 0;
        my $Da_value = 0;
        my $Cp_value = 0;
        my $Neis_D_value = 0;
        
        
        if ($Adenine_count_fixed <= $minimum_allele_count){
            $Adenine_count_fixed = 0
        }
        if ($Thymine_count_fixed <= $minimum_allele_count){
            $Thymine_count_fixed = 0
        }
        if ($Cytosine_count_fixed <= $minimum_allele_count){
            $Cytosine_count_fixed = 0
        }
        if ($Guanine_count_fixed <= $minimum_allele_count){
            $Guanine_count_fixed = 0
        }
        if ($Deletion_count_fixed <= $minimum_allele_count){
            $Deletion_count_fixed = 0
        }
        if ($Adenine_count_polymorphic <= $minimum_allele_count){
            $Adenine_count_polymorphic = 0
        }
        if ($Thymine_count_polymorphic <= $minimum_allele_count){
            $Thymine_count_polymorphic = 0
        }
        if ($Cytosine_count_polymorphic <= $minimum_allele_count){
            $Cytosine_count_polymorphic = 0
        }
        if ($Guanine_count_polymorphic <= $minimum_allele_count){
            $Guanine_count_polymorphic = 0
        }
        if ($Deletion_count_polymorphic <= $minimum_allele_count){
            $Deletion_count_polymorphic = 0
        }
        
        my $read_nucleotides_fixed = $Adenine_count_fixed + $Thymine_count_fixed + $Cytosine_count_fixed + $Guanine_count_fixed + $Deletion_count_fixed;
        my $read_nucleotides_polymorphic = $Adenine_count_polymorphic + $Thymine_count_polymorphic + $Cytosine_count_polymorphic + $Guanine_count_polymorphic + $Deletion_count_polymorphic;
        
        if ($read_nucleotides_fixed > $maximum_read_depth || $read_nucleotides_polymorphic > $maximum_read_depth){
            next;
        }
        
        if ($read_nucleotides_fixed < $minimum_read_depth || $read_nucleotides_polymorphic < $minimum_read_depth){
            next;
        }
        
        if ($start_position_boolean eq "false"){
            $window_first_position = $array_of_line[1];
            $start_position_boolean = "true"
        }
        
        if ($Scaffold ne "$old_scaffold"){
            if ($window_counter!= 0 && $window_counter >= 0.1 * $sex_SNP_finder_window_size){
                $snp_counter = $snp_counter / $window_counter * $sex_SNP_finder_window_size;
                print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n"
            }
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
        
        $old_scaffold = $Scaffold;
        $window_last_position = $array_of_line[1];
        
        
        if ($read_nucleotides_fixed >= $minimum_read_depth){
            if($read_nucleotides_polymorphic >= $minimum_read_depth){
                # window counter
                $window_counter = $window_counter + 1;
                
                # Adenine
                
                my $Adenine_proportion_fixed = $Adenine_count_fixed / $read_nucleotides_fixed;
                my $Non_adenine_proportion_polymorphic = 1 - $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_adenine_proportion_polymorphic > 0 && $Adenine_proportion_fixed <= 1){
                    if($Non_adenine_proportion_polymorphic <= 1 && $Adenine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Adenine_proportion_fixed >= $fixed_threshold){
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                        }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                        }
                }
                
                
                #Thymine
                
                my $Thymine_proportion_fixed = $Thymine_count_fixed / $read_nucleotides_fixed;
                my $Non_thymine_proportion_polymorphic = 1 - $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_thymine_proportion_polymorphic > 0 && $Thymine_proportion_fixed <= 1){
                    if($Non_thymine_proportion_polymorphic <= 1 && $Thymine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Thymine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                        }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                        }
                }
                
                # Cytosine
                
                my $Cytosine_proportion_fixed = $Cytosine_count_fixed / $read_nucleotides_fixed;
                my $Non_cytosine_proportion_polymorphic = 1 - $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_cytosine_proportion_polymorphic > 0 && $Cytosine_proportion_fixed <= 1){
                    if($Non_cytosine_proportion_polymorphic <= 1 && $Cytosine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Cytosine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                        }
                }
                
                # Guanine
                
                my $Guanine_proportion_fixed = $Guanine_count_fixed / $read_nucleotides_fixed;
                my $Non_guanine_proportion_polymorphic = 1 - $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_guanine_proportion_polymorphic > 0 && $Guanine_proportion_fixed <= 1 ){
                    if ($Non_guanine_proportion_polymorphic <= 1 && $Guanine_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Guanine_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Deletion_proportion_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                    }
                    if ($Deletion_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Deletion_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Deletion_proportion_polymorphic\n";
                        }
                }
                
                # Deletion
                
                my $Deletion_proportion_fixed = $Deletion_count_fixed / $read_nucleotides_fixed;
                my $Non_deletion_proportion_polymorphic = 1 - $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
                if ($Non_deletion_proportion_polymorphic > 0 && $Deletion_proportion_fixed <= 1){
                    if ($Non_deletion_proportion_polymorphic <= 1 && $Deletion_proportion_fixed > 0){
                        ($Fst_value, $Pi_fixed, $Pi_polymorphic) = Fst($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
                        $Fst_boolean = $Fst_boolean + 1;
                        ($Dxy_value, $Cp_value, $Neis_D_value) = Dxy($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $maximum_neis_D);
                        $Dxy_boolean = $Dxy_boolean + 1;
                        $Da_value = $Dxy_value - ($Pi_fixed + $Pi_polymorphic) / 2;
                    }
                }
                
                if ($Deletion_proportion_fixed >= $fixed_threshold){
                    my $Adenine_proportion_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Thymine_proportion_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Cytosine_proportion_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
                    my $Guanine_proportion_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
                    if ($Adenine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Adenine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Adenine_proportion_polymorphic\n";
                        }
                    if ($Thymine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Thymine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Thymine_proportion_polymorphic\n";
                        }
                    if ($Cytosine_proportion_polymorphic <= $maximum_polymorphic_frequency && $Cytosine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                        # snp counter
                        $snp_counter = $snp_counter + 1;
                        print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cytosine_proportion_polymorphic\n";
                    }
                    if ($Guanine_proportion_polymorphic <= $maximum_polymorphic_frequency &&
                        $Guanine_proportion_polymorphic >= $minimum_polymorphic_frequency){
                            # snp counter
                            $snp_counter = $snp_counter + 1;
                            print $OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Guanine_proportion_polymorphic\n";
                        }
                }
            }
        }
        print $FST_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Fst_value\n";
        print $DXY_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Dxy_value\n";
        print $DA_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Da_value\n";
        print $CP_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Cp_value\n";
        print $NEIS_D_OUTPUT "$Scaffold\t$Start_position\t$End_position\t$Feature\t$Neis_D_value\n";
        $i = $i + 1;
        my $d = 1000000;
        if ($i % $d == 0){
            print "[$i] nucleotide positions analyzed.\n";
        }
        if ($window_counter == $sex_SNP_finder_window_size){
            print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n";
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
    }
}

else {
    print "ERROR: You need to declare which pool is fixed. Enter either pool1 or pool2 into the --fixed_population option.\n"
}

close $INPUT;
close $OUTPUT;
close $SL_OUTPUT;

sub Pi
{
    my ($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed, $Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic) = @_;
    my ($Adenine_frequency_fixed, $Thymine_frequency_fixed, $Cytosine_frequency_fixed, $Guanine_frequency_fixed, $Deletion_frequency_fixed, $Adenine_frequency_polymorphic, $Thymine_frequency_polymorphic, $Cytosine_frequency_polymorphic, $Guanine_frequency_polymorphic, $Deletion_frequency_polymorphic, $Adenine_frequency_average, $Thymine_frequency_average, $Cytosine_frequency_average, $Guanine_frequency_average, $Deletion_frequency_average) = Frequencies($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
    # Pi calculations
    my $Pi_fixed = (1 - $Adenine_frequency_fixed ** 2 - $Thymine_frequency_fixed ** 2 - $Cytosine_frequency_fixed ** 2 - $Guanine_frequency_fixed ** 2 - $Deletion_frequency_fixed ** 2);
    my $Pi_polymorphic = (1 - $Adenine_frequency_polymorphic ** 2 - $Thymine_frequency_polymorphic ** 2 - $Cytosine_frequency_polymorphic ** 2 - $Guanine_frequency_polymorphic ** 2 - $Deletion_frequency_polymorphic ** 2);
    my $Pi_total = (1 - $Adenine_frequency_average ** 2 - $Thymine_frequency_average ** 2 - $Cytosine_frequency_average ** 2 - $Guanine_frequency_average ** 2 - $Deletion_frequency_average ** 2);
    return($Pi_fixed, $Pi_polymorphic, $Pi_total);
    exit;
}

sub Frequencies
{
    my ($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic) = @_;
    # Fixed frequencies
    my $Adenine_frequency_fixed = $Adenine_count_fixed / $read_nucleotides_fixed;
    my $Thymine_frequency_fixed = $Thymine_count_fixed / $read_nucleotides_fixed;
    my $Cytosine_frequency_fixed = $Cytosine_count_fixed / $read_nucleotides_fixed;
    my $Guanine_frequency_fixed = $Guanine_count_fixed / $read_nucleotides_fixed;
    my $Deletion_frequency_fixed = $Deletion_count_fixed / $read_nucleotides_fixed;
    # Polymorphic frequencies
    my $Adenine_frequency_polymorphic = $Adenine_count_polymorphic / $read_nucleotides_polymorphic;
    my $Thymine_frequency_polymorphic = $Thymine_count_polymorphic / $read_nucleotides_polymorphic;
    my $Cytosine_frequency_polymorphic = $Cytosine_count_polymorphic / $read_nucleotides_polymorphic;
    my $Guanine_frequency_polymorphic = $Guanine_count_polymorphic / $read_nucleotides_polymorphic;
    my $Deletion_frequency_polymorphic = $Deletion_count_polymorphic / $read_nucleotides_polymorphic;
    # Average frequencies
    my $Adenine_frequency_average = ($Adenine_frequency_fixed + $Adenine_frequency_polymorphic) / 2;
    my $Thymine_frequency_average = ($Thymine_frequency_fixed + $Thymine_frequency_polymorphic) / 2;
    my $Cytosine_frequency_average = ($Cytosine_frequency_fixed + $Cytosine_frequency_polymorphic) / 2;
    my $Guanine_frequency_average = ($Guanine_frequency_fixed + $Guanine_frequency_polymorphic) / 2;
    my $Deletion_frequency_average = ($Deletion_frequency_fixed + $Deletion_frequency_polymorphic) / 2;
    return ($Adenine_frequency_fixed, $Thymine_frequency_fixed, $Cytosine_frequency_fixed, $Guanine_frequency_fixed, $Deletion_frequency_fixed, $Adenine_frequency_polymorphic, $Thymine_frequency_polymorphic, $Cytosine_frequency_polymorphic, $Guanine_frequency_polymorphic, $Deletion_frequency_polymorphic, $Adenine_frequency_average, $Thymine_frequency_average, $Cytosine_frequency_average, $Guanine_frequency_average, $Deletion_frequency_average);
}

sub Fst
{
    my ($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic) = @_;
    my ($Pi_fixed, $Pi_polymorphic, $Pi_total) = Pi($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
    my $Pi_fixed_fst = $Pi_fixed * ($read_nucleotides_fixed / ($read_nucleotides_fixed - 1));
    my $Pi_polymorphic_fst = $Pi_polymorphic * ($read_nucleotides_polymorphic / ($read_nucleotides_polymorphic - 1));
    my $smaller_read_count = min($read_nucleotides_fixed, $read_nucleotides_polymorphic);
    my $Pi_total_fst = $Pi_total * ($smaller_read_count / ($smaller_read_count - 1));
    my $Pi_within = ($Pi_fixed_fst + $Pi_polymorphic_fst) / 2;
    my $Fst = ($Pi_total_fst - $Pi_within) / $Pi_total_fst;
    return ($Fst, $Pi_fixed, $Pi_polymorphic);
    exit;
}

sub Dxy
{
    my ($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed, $Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic, $max_neis_D) = @_;
    my ($Adenine_frequency_fixed, $Thymine_frequency_fixed, $Cytosine_frequency_fixed, $Guanine_frequency_fixed, $Deletion_frequency_fixed, $Adenine_frequency_polymorphic, $Thymine_frequency_polymorphic, $Cytosine_frequency_polymorphic, $Guanine_frequency_polymorphic, $Deletion_frequency_polymorphic, $Adenine_frequency_average, $Thymine_frequency_average, $Cytosine_frequency_average, $Guanine_frequency_average, $Deletion_frequency_average) = Frequencies($Adenine_count_fixed, $Thymine_count_fixed, $Cytosine_count_fixed,$Guanine_count_fixed, $Deletion_count_fixed, $read_nucleotides_fixed, $Adenine_count_polymorphic, $Thymine_count_polymorphic, $Cytosine_count_polymorphic, $Guanine_count_polymorphic, $Deletion_count_polymorphic, $read_nucleotides_polymorphic);
    my $Dxy = ($Adenine_frequency_fixed * $Thymine_frequency_polymorphic) + ($Adenine_frequency_fixed * $Cytosine_frequency_polymorphic) + ($Adenine_frequency_fixed * $Guanine_frequency_polymorphic) + ($Adenine_frequency_fixed * $Deletion_frequency_polymorphic) + ($Thymine_frequency_fixed * $Adenine_frequency_polymorphic) + ($Thymine_frequency_fixed * $Cytosine_frequency_polymorphic) + ($Thymine_frequency_fixed * $Guanine_frequency_polymorphic) + ($Thymine_frequency_fixed * $Deletion_frequency_polymorphic) + ($Cytosine_frequency_fixed * $Adenine_frequency_polymorphic) + ($Cytosine_frequency_fixed * $Thymine_frequency_polymorphic) + ($Cytosine_frequency_fixed * $Guanine_frequency_polymorphic) + ($Cytosine_frequency_fixed * $Deletion_frequency_polymorphic) + ($Guanine_frequency_fixed * $Adenine_frequency_polymorphic) + ($Guanine_frequency_fixed * $Thymine_frequency_polymorphic) + ($Guanine_frequency_fixed * $Cytosine_frequency_polymorphic) + ($Guanine_frequency_fixed * $Deletion_frequency_polymorphic) + ($Deletion_frequency_fixed * $Adenine_frequency_polymorphic) + ($Deletion_frequency_fixed * $Thymine_frequency_polymorphic) + ($Deletion_frequency_fixed * $Cytosine_frequency_polymorphic) + ($Deletion_frequency_fixed * $Guanine_frequency_polymorphic);
    my $Cp = 0.5 * (abs($Adenine_frequency_fixed - $Adenine_frequency_polymorphic) + abs($Thymine_frequency_fixed - $Thymine_frequency_polymorphic) + abs($Cytosine_frequency_fixed - $Cytosine_frequency_polymorphic) + abs($Guanine_frequency_fixed - $Guanine_frequency_polymorphic) + abs($Deletion_frequency_fixed - $Deletion_frequency_polymorphic));
    my $Jxy = $Adenine_frequency_fixed * $Adenine_frequency_polymorphic + $Thymine_frequency_fixed * $Thymine_frequency_polymorphic + $Cytosine_frequency_fixed * $Cytosine_frequency_polymorphic + $Guanine_frequency_fixed * $Guanine_frequency_polymorphic + $Deletion_frequency_fixed * $Deletion_frequency_polymorphic;
    my $Jx = $Adenine_frequency_fixed ** 2 + $Thymine_frequency_fixed ** 2 + $Cytosine_frequency_fixed ** 2 + $Guanine_frequency_fixed ** 2 + $Deletion_frequency_fixed ** 2;
    my $Jy = $Adenine_frequency_polymorphic ** 2 + $Thymine_frequency_polymorphic ** 2 + $Cytosine_frequency_polymorphic ** 2 + $Guanine_frequency_polymorphic ** 2 + $Deletion_frequency_polymorphic ** 2;
    my $I = $Jxy / (($Jx * $Jy) ** 0.5 );
    my $neis_D = 0;
    if ($I == 0){
        $neis_D = $max_neis_D;
    }
    else {
        $neis_D = -1 * ln($I);
    }
    return ($Dxy, $Cp, $neis_D);
}

sub Usage
{
    my $command = $0;
    $command =~ s#^[^\s]/##;
    printf STDERR "@_\n" if ( @_ );
    printf STDERR "\nperl Sex_SNP_finder_GA.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --minimum_read_depth=[value greater than 0] --maximum_read_depth=[value greater than the minimum_read_depth] --minimum_allele_count=[minimum number of allele counts at a site before an allele is considered real and not a sequencing error] --sex_SNP_finder_window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --fst_output_file=fst_output_file.igv --dxy_output_file=dxy_output_file.igv --da_output_file=da_output_file.igv --cp_output_file=cp_output_file.igv --neis_D_output_file=neis_D_output_file.igv--description=description of file to be used in IGV header [--help|-?]\n\n";
    exit;
}
