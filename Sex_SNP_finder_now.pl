#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=pod
 
=head1 NAME
 
 Sex_SNP_finder_now.pl
 
=head1 AUTHORS
 
 Will Gammerdinger - Program Designer and Matt Conte - Coding Advisor
 
=head1 DESCRIPTION
 
 Sex_SNP_finder_now.pl was a derivation from Sex_SNP_finder.pl. It is designed to look for regions of differentiation using a non overlapping window approach. Sex_SNP_finder_now.pl was designed to find nucleotides that were fixed at a site in one population and polymorphic at the same site in a different population. Initially, the program was used to compare between males and females, in order to find mutations responsible for sex-determination. However, it is applicable to any population where one nucleotide is fixed in one population and polymorphic in the other.
 
=head1 EXAMPLE
 
 The format should be perl Sex_SNP_finder_now.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --read_depth=[value greater than 0] --window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --description=description of file to be used in IGV header [--help|-?]\n\n
 
=head1 LICENSE
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=head1 VERSION
 
 version0.0.6 - The biallelic assumption was corrected. Additionally, the non-overlapping window output was changed to an IGV readable output.
 
 version0.0.5 - Names of Sex_finder software were changed from Sex_finder to Sex_SNP_finder. Added a non-overlapping window approach to Sex_SNP_finder.pl to make Sex_SNP_finder_now.pl. This non-overlapping window skips sites that fail to meet the read depth threshold and gaps.
 
 version0.0.4 - Returned a guide prompt when option fields were left empty.
 
 version0.0.3 - Changed the input and output from arguements to options.
 
 version0.0.2 - Allowed declaration of pooled data set
 
 version0.0.1 - Original Version
 
=cut


my ($input_file, $output_file, $fixed_population, $fixed_threshold, $minimum_polymorphic_frequency, $maximum_polymorphic_frequency, $read_depth, $window_size, $non_overlapping_window_output_file, $description, $help) = ("empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty", "empty");

GetOptions(
"input_file=s"                          => \$input_file,
"output_file=s"                         => \$output_file,
"fixed_population=s"                    => \$fixed_population,
"fixed_threshold=s"                     => \$fixed_threshold,
"minimum_polymorphic_frequency=s"       => \$minimum_polymorphic_frequency,
"maximum_polymorphic_frequency=s"       => \$maximum_polymorphic_frequency,
"read_depth=s"                          => \$read_depth,
"window_size=s"                         => \$window_size,
"non_overlapping_window_output_file=s"  => \$non_overlapping_window_output_file,
"description=s"                         => \$description,
"help|?"                                => \$help
) or Usage ( "Invalid command-line option.");

Usage() if defined $help;

if ($input_file eq "empty" || $output_file eq "empty" || $fixed_population eq "empty" || $fixed_threshold eq "empty" || $minimum_polymorphic_frequency eq "empty" || $maximum_polymorphic_frequency eq "empty" || $read_depth eq "empty" || $window_size eq "empty" || $non_overlapping_window_output_file eq "empty" || $description eq "empty"){
    die "\nERROR: The format should be perl Sex_SNP_finder_now.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --read_depth=[value greater than 0] --window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --description=description of file to be used in IGV header [--help|-?]\n\nOne or more of your option fields is empty.\n\nFor more information, use the command perldoc Sex_SNP_finder_now.pl\n\n"
}

open (my $INPUT, "<$input_file");

open (my $OUTPUT, ">$output_file");
open (my $SL_OUTPUT, ">$non_overlapping_window_output_file");


print "Your fixed population is ", $fixed_population, ".\n";
print "Your minimum read depth is ", $read_depth, " nucleotides.\n";
print "Your minimum fixed threshold is ", $fixed_threshold, ".\n";
print "Your minumum polymorphic frequency is ", $minimum_polymorphic_frequency, ".\n";
print "Your maximum polymorphic frequency is ", $maximum_polymorphic_frequency, ".\n";

if ($minimum_polymorphic_frequency > $maximum_polymorphic_frequency){
    die "ERROR: The minimum polymorphic frequency is greater than the maximum polymorphic frequency.\nRemember that these allele frequencies of all sites that are not the same as the fixed allele.\n"
}

print "Your window size is ", $window_size, ".\n";


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

my $i = 0;
my $window_counter=0;
my $snp_counter=0;
my $old_scaffold = $first_linkage_group;
my $start_position_boolean = "false";
my $window_first_position;
my $window_last_position;

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
        
        my $read_nucleotides_fixed = $Adenine_count_fixed + $Thymine_count_fixed + $Cytosine_count_fixed + $Guanine_count_fixed + $Deletion_count_fixed;
        my $read_nucleotides_polymorphic = $Adenine_count_polymorphic + $Thymine_count_polymorphic + $Cytosine_count_polymorphic + $Guanine_count_polymorphic + $Deletion_count_polymorphic;
        
        if ($start_position_boolean eq "false"){
            $window_first_position = $array_of_line[1];
            $start_position_boolean = "true"
        }
        
        if ($Scaffold ne "$old_scaffold"){
            if ($window_counter!= 0 && $window_counter >= 0.1 * $window_size){
                $snp_counter = $snp_counter / $window_counter * $window_size;
                print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n"
            }
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
        
        $old_scaffold = $Scaffold;
        $window_last_position = $array_of_line[1];
        
        
        if ($read_nucleotides_fixed >= $read_depth){
            if($read_nucleotides_polymorphic >= $read_depth){
                # window counter
                $window_counter = $window_counter + 1;

                # Adenine
                
                my $Adenine_proportion_fixed = $Adenine_count_fixed / $read_nucleotides_fixed;
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
        $i = $i + 1;
        my $d = 1000000;
        if ($i % $d == 0){
            print "[$i] nucleotide positions analyzed.\n";
        }

        if ($window_counter == $window_size){
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
        
        my $read_nucleotides_fixed = $Adenine_count_fixed + $Thymine_count_fixed + $Cytosine_count_fixed + $Guanine_count_fixed + $Deletion_count_fixed;
        my $read_nucleotides_polymorphic = $Adenine_count_polymorphic + $Thymine_count_polymorphic + $Cytosine_count_polymorphic + $Guanine_count_polymorphic + $Deletion_count_polymorphic;
        
        if ($start_position_boolean eq "false"){
            $window_first_position = $array_of_line[1];
            $start_position_boolean = "true"
        }
        
        if ($Scaffold ne "$old_scaffold"){
            if ($window_counter!= 0 && $window_counter >= 0.1 * $window_size){
                $snp_counter = $snp_counter / $window_counter * $window_size;
                print $SL_OUTPUT "$Scaffold\t$window_first_position\t$window_last_position\t$SL_feature\t$snp_counter\n"
            }
            $window_counter = 0;
            $snp_counter = 0;
            $start_position_boolean = "false";
        }
        
        $old_scaffold = $Scaffold;
        $window_last_position = $array_of_line[1];

        
        if ($read_nucleotides_fixed >= $read_depth){
            if($read_nucleotides_polymorphic >= $read_depth){
                
                # window counter
                $window_counter = $window_counter + 1;
                
                # Adenine
                
                my $Adenine_proportion_fixed = $Adenine_count_fixed / $read_nucleotides_fixed;
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
        $i = $i + 1;
        my $d = 1000000;
        if ($i % $d == 0){
            print "[$i] nucleotide positions analyzed.\n";
        }
        
        if ($window_counter == $window_size){
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

sub Usage
{
    my $command = $0;
    $command =~ s#^[^\s]/##;
    printf STDERR "@_\n" if ( @_ );
    printf STDERR "\nThe format should be perl Sex_SNP_finder_now.pl --input_file=input_file.sync --output_file=output_file.igv --fixed_population=pool[1 or 2] --fixed_threshold=[value between 0 and 1] --minimum_polymorphic_frequency=[value between 0 and 1] --maximum_polymorphic_frequency=[value between 0 and 1] --read_depth=[value greater than 0] --window_size=[value greater than 0] --non_overlapping_window_output_file=non_overlapping_window_output_file.igv --description=description of file to be used in IGV header [--help|-?]\n\nOne or more of your option fields is empty.\n\nFor more information, use the command perldoc Sex_SNP_finder_now.pl\n\n";
    exit;
}
