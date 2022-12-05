#_____________________________________________CONSERVATION________________________________________________________________________________



use warnings;
use strict;
#use diagnostics;

use File::Path qw(make_path);

use List::Util qw(sum uniq any max reduce);
use List::Compare;

use JSON;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );

my $info = Sys::Info->new;
my $cpu  = $info->device( 'CPU' );
my $NumberCPU = $cpu->count;
my $NumberCPUx2 = $NumberCPU * 2;

use Parallel::ForkManager;
my $forks = $NumberCPUx2 or die "Usage: $0 N\n";

use Storable qw(dclone);





########### COMPUTATION OF THE CONSERVATION ACROSS DIFFERENT LEVEL OF CLASSIFICATION FOR EACH CALCULTED REAL INTERACTION ###########



my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



#my @Flu = ("FluX");
#my @StructureStates = ("Prefusion");
    my @Flu = ("FluA", "FluB");
    my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");

    my @AllInteractionTypes = ("All", "Ionic", "Repulsive", "RepulsivePositive", "RepulsiveNegative", "Hydrogen", "Hydrogen_SD_SD", "Hydrogen_SD_C", "Hydrophobic");

    my @Indications = ( "Total", "Intra", "Inter", "HA1-HA1", "HA1-HA2", "HA2-HA2",
                          "IntraHA1-HA1", "IntraHA1-HA2", "IntraHA2-HA2",
                          "InterHA1-HA1", "InterHA1-HA2", "InterHA2-HA2" );

    my @InteractionTypes = qw(Inter-Ionic             Intra-Ionic

                           Inter-RepulsivePositive    Intra-RepulsivePositive
                           Inter-RepulsiveNegative    Intra-RepulsiveNegative

                           Inter-Hydrogen_SD_SD       Intra-Hydrogen_SD_SD
                           Inter-Hydrogen_SD_C        Intra-Hydrogen_SD_C
                           Inter-Hydrophobic          Intra-Hydrophobic);



   print "...\nWill process all following interaction types :\n";
   print "$_\n" foreach (@AllInteractionTypes);
   print "\n.........\n !!!! Change the array named AllInteractionsTypes if we add new interaction types !!!! \n.........\n\n";
   print "\n.........\nMachine CPUs = $NumberCPU\n.........\n\n";
   print "\n.........\nCheck subs interact_symetry_check get_all_corresponding_interaction_and_assess_similarity & keep_one_interaction if wanna assess interaction assymetry\n.........\n";


















############### GET HASH WITH ALL FLU / STRUCTURE STATE / SUBTYPE AND CODE PDB ###############
      my %H_state_subtype_code = ();
      my @List_all_codePDB = ();

            foreach my $flu (@Flu) {
              foreach my $StructureState (@StructureStates){
               unless ($flu eq "FluB" and $StructureState eq "Intermediates" or $flu eq "FluB" and $StructureState eq "Postfusion"){

                my $dirInputSubtype = "$path_environment/Interactions/Result/$flu/$StructureState";
                my @SubtypesInputToProcess = opendir_readdir ($dirInputSubtype);

                foreach my $Subtype (@SubtypesInputToProcess) {

                  next if ($Subtype =~ m/^\./);
                  my $PathToCodesPDB = $dirInputSubtype . "/" . $Subtype;
                  my @CodesPDBInputToProcess = opendir_readdir ($PathToCodesPDB);

                  foreach my $CodePDB (@CodesPDBInputToProcess) {
                    next if ($CodePDB =~ m/^\./);
                    push @List_all_codePDB, $CodePDB;
                    push(@{ $H_state_subtype_code{$flu}{$StructureState}{$Subtype} }, $CodePDB);

                    }
                  }
                 }
                }
              }

      print "\n\n\n";
      print Dumper \%H_state_subtype_code;
      print "\n\n\nABOVE : Hash despicting all structures for conservation calculation\n\n\n";
      print "\n\n\n" x 2;
      my $total_codes_pdb = scalar(@List_all_codePDB);
############### GET HASH WITH ALL FLU / STRUCTURE STATE / SUBTYPE AND CODE PDB ###############






















 ############### GET HASH ALL UNIQUE INTERACTIONS (NON DUPLICATED) FOR EACH FLU/STATE SECTION : UNIQUE INTERACTION ###############
 my %Full_H_all_interactions = ();
 my $count_pdb = 0;

 print uc("Screening uniques interactions amongst structures...\n");
 print uc("--> total PDB codes to parse : $total_codes_pdb\n\n\n");

 foreach my $flu (@Flu) {
   print "$flu\n\n";


   foreach my $State (@StructureStates) {
    unless ($flu eq "FluB" and $State eq "Intermediates" or $flu eq "FluB" and $State eq "Postfusion"){
      print "$State\n";
      my %H_all_interactions = ();


     foreach my $Subtype ( keys %{ $H_state_subtype_code{$flu}{$State} } ) {
       print "\n\n$Subtype\n";


       foreach my $CodePDB ( @{ $H_state_subtype_code{$flu}{$State}{$Subtype} } ) {
         $count_pdb++;
         print "  $count_pdb / $total_codes_pdb  ";


         foreach my $Type_interaction (@InteractionTypes) {

           my $Path = "$path_environment/Interactions/Result/$flu/$State/$Subtype/$CodePDB/";
           my $file = $Type_interaction . "_" . $CodePDB . ".txt";
           open (INPUTT   , "<", "$Path$file") or die "Pb d'ouverture $Path$file : $!";


           while (<INPUTT>) {
              # H1N1 1RU7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299
              #  0    1      2        3   4   5  6    7       8        9     10   11 12  13  14     15     16      17

              chomp;

              my @Splt = split(" ", $_);
              foreach (@Splt) { $_ =~ s/^\s+|\s+$//g; }

              my $potential_unique_order   = join(" ", $Splt[2], $Splt[4], $Splt[5], $Splt[6], $Splt[11], $Splt[12], $Splt[13]);
              my $potential_unique_reverse = join(" ", $Splt[2], $Splt[11], $Splt[12], $Splt[13], $Splt[4], $Splt[5], $Splt[6]);


              #CHECK IF INTERACTION (OR REVERSE ORDER) IS ALREADY PRESENT IN OUR HASH
              unless ( grep{ $_ =~ /$potential_unique_order|$potential_unique_reverse/ } @{ $H_all_interactions{$Type_interaction} }  )   {

                      push @{ $H_all_interactions{$Type_interaction} }, $potential_unique_order;

                                    }
                                  }
                close INPUTT;
                    }
                  }
                }

                my $ref_h = \%H_all_interactions;
                $Full_H_all_interactions{$flu}{$State} = $ref_h;

             }
           }
         }

 print "\n\n\n\n";
 print " ---> All unique interactions amongst all structures have been retrieved.\n\n\n\n";

 my $path_uniques = "$path_environment/Interactions/Result/Conservation";

 my $Path_json_all_uniques = h_json_and_dumper (\%Full_H_all_interactions, $path_uniques, "All_Unique_Interactions");

 undef %Full_H_all_interactions;
 ############### GET HASH ALL UNIQUE INTERACTIONS (NON DUPLICATED) FOR EACH FLU/STATE SECTION : UNIQUE INTERACTION ###############













  ############### CONSERVATION PART ###############
  ############### CONSERVATION PART ###############
  ############### CONSERVATION PART ###############
  print "\n\n" x 5;
  print uc("__ Start of conservation calculation...\n\n");



  #WORK WITH HASH IN HARDWARE
  my $href_all_uniques = get_json_data_hardsaved ($Path_json_all_uniques);
  my @all_unique_conserved = ();



 foreach my $flu (@Flu) {

  foreach my $State (@StructureStates) {
   unless ($flu eq "FluB" and $State eq "Intermediates" or $flu eq "FluB" and $State eq "Postfusion"){

    #key 1 : type interaction, key 2 : Unique interaction, value = array ref with all correspondings
    my $ref_H_uniq_interact = $href_all_uniques->{$flu}{$State};

    #FOREACH UNIQUE INTERACTION WE ADD ALL THE CONRRESPONDING INTERACTIONS FOUND
    my %All_corresponding_interactions = ();


          foreach my $Type_interaction (@InteractionTypes) {

            if ( exists($ref_H_uniq_interact->{$Type_interaction}) ) {

              #New Hash foreach interaction type
              my %H_conservation = ();


              print "\n\n\n...............................................................................................................................\n";
              my $total_uniqs = scalar (@{ $ref_H_uniq_interact->{$Type_interaction} });
              my $count_uniq = 0;
              print "\n\n$flu $State $Type_interaction : nb uniqs to proceed (get_all_corresponding_interaction_and_assess_similarity) = $total_uniqs\n\n";



              #PARALLEL PROCESSING OF CONSERVATION ASSESSMENTS
              my $pm = Parallel::ForkManager->new($NumberCPU);

              $pm -> run_on_finish ( # called BEFORE the first call to start()
                sub {
                  my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

                  # retrieve data structure from child
                  if (defined($data_structure_reference)) {  # children are not forced to send anything


                  foreach my $type (keys %{ $data_structure_reference } ) {
                    foreach my $uniq (keys %{ $data_structure_reference->{$type} } ) {
                      foreach my $pushed ( @{ $data_structure_reference->{$type}{$uniq} } ) {

                         push @{ $All_corresponding_interactions{$type}{$uniq} }, $pushed;

                        }
                      }
                    }

                  }
                }
              );

              $pm->run_on_wait( sub {
                print " "
              },
              0.5
              );


              #START PARALLELIZATION to grab all interaction CORRESPONDING TO EACH UNIQUE INTERACTIONS
              LOOP:foreach my $Unique_interaction (@{ $ref_H_uniq_interact->{$Type_interaction} }) {

                      $count_uniq++;
                      print "$count_uniq ";

                      my $pid = $pm->start and next LOOP;

                      #Get all interaction corresponding to a specific unique interaction (More comments within function)
                      #Find all interaction corresponding to "Unique_interaction" and assess similarity level to ADD CONSERVATION INFOS
                      my $ref_h = get_all_corresponding_interaction_and_assess_similarity ($flu, $State, $Type_interaction, $Unique_interaction, \%H_state_subtype_code);

                      $pm->finish(0, $ref_h);
                    }

              print "\n\n";
              $pm->wait_all_children;



               my $OutP = "$path_environment/Interactions/Result/Conservation/$flu/$State/$Type_interaction";
               my $Path_json_all_corresponding = h_json_and_dumper (\%All_corresponding_interactions, $OutP, "AllInteraction_and_correspondance_$Type_interaction");

               undef %All_corresponding_interactions;
               my $h_all_corresponding_interacts = get_json_data_hardsaved ($Path_json_all_corresponding);



               print "\n\n ---> All interactions corresponding to each unique have been retrieved for every --- $flu $State $Type_interaction --- structures\n";
               print "      Conservation info have been added\n\n";





               my %h_conserved_interacts = ();
               # IF WANNA ASSES ASSYMETRY : LOOK HERE
               #For each interaction, amongst all corresponding (conserved) interaction --> keep THE MOST RELEVANT per codepdb

               my $nb_uniq_h_all_corresponding_interact = scalar(keys %{ $h_all_corresponding_interacts->{$Type_interaction} });
               my $nb_proceed = 0;
               print "\n\n$flu $State $Type_interaction : nb uniqs to proceed (keep_mostsimilar_interaction) = $nb_uniq_h_all_corresponding_interact\n\n";

               foreach my $Uniq_interaction (keys %{ $h_all_corresponding_interacts->{$Type_interaction} }) {

                  $nb_proceed++;
                  print " $nb_proceed";

                  my @definitive_keeps = ();

                  #SEGREGATE CORRESPONDING INTERACTIONS DEPENDING ON THE ADDED CONSERVATION INFOS
                  foreach my $PDBcode (@List_all_codePDB) {

                      my @Grep = grep{ /$PDBcode/ } @{ $h_all_corresponding_interacts->{$Type_interaction}{$Uniq_interaction} };

                      my @Grep_1 = grep{ /Identical_Strictly/ }  @Grep;        my @Grep_2 = grep{ /TypeConserved_Strictly/ }  @Grep;
                      my @Grep_3 = grep{ /Identical_InTheArea/ } @Grep;        my @Grep_4 = grep{ /TypeConserved_InTheArea/ } @Grep;


                      my $definitive_keep;
                      #We keep the most similar interaction to $Uniq_interaction and assess for symmetry
                      WLOOP:
                      while (@Grep) {

                          my ($Keep, $ref_grep_keep) = keep_mostsimilar_interaction (\@Grep_1, \@Grep_2, \@Grep_3, \@Grep_4, $Uniq_interaction);
                          #TypeConserved_InTheArea H18N11	4mc5	Inter-Hydrophobic	CE2 TYR 458 HA2 Proto1 old=423   |     CD1 LEU 463 HA2 Proto3 old=428 ->   3.96914474414828
                          #             0            1      2            3          4  5   6   7   8      9        10    11   12  13  14   15      16    17       18

                          my $defo_keep = interact_symmetry_check ($Keep, $ref_grep_keep, \@Grep);
                          $definitive_keep = $defo_keep if defined($defo_keep);

                          last WLOOP if defined($definitive_keep);

                        }

                      push @definitive_keeps, $definitive_keep if defined($definitive_keep);

                 }

                  $h_conserved_interacts{$Type_interaction}{$Uniq_interaction} = \@definitive_keeps if @definitive_keeps;

               }

               my $OutConserv = "$path_environment/Interactions/Result/Conservation/$flu/$State/$Type_interaction";
               print "\n\n";
               my $path_h_conserved_interacts = h_json_and_dumper (\%h_conserved_interacts, $OutConserv, "Interaction_Conserv_$Type_interaction");

               print "\n\n ---> Conserved symetric interactions have been retrieved for every --- $flu $State $Type_interaction --- structures\n";

               undef %h_conserved_interacts;
               undef $h_all_corresponding_interacts;
               print "\n\n" x 3;
               print uc("\n\n\n___ Conserved ** $Type_interaction ** interactions corresponding to each unique one have been selected for each code PDB ___\n\n\n");

               print "\n\n\n\n";






           #ASSESS CONSERVATION OF SPECIFIC INTERACTION ACCROSS CLASSIFICATIONS
           if (-e $path_h_conserved_interacts) {

             my $ref_h_all_corresp = get_json_data_hardsaved ($path_h_conserved_interacts);


             my $State_data_array = "Prefusion";    $State_data_array = "Postfusion_Intermediary" if $State =~ m/Intermediates|Postfusion/;
             my $Path_class = "$path_environment/Interactions/Data/ArrayTableauPDB/$flu/$State_data_array/HashAllConditionsGrep_$State_data_array.dat";

             my $ref_H_category = get_json_data_hardsaved ($Path_class);
             delete $ref_H_category->{"Subtypes_Hosts"} if exists $ref_H_category->{"Subtypes_Hosts"};

             next if !defined($ref_H_category);
             my $copy_H_category = dclone($ref_H_category);

             print "\n\nAssessing conservation for different classifications :\n\n";


             while (keys(%$copy_H_category)) {

             #Recursive function to get all groups conservation IE CladeH1
               my ($last_key, $ref_array_v, $ref_key_list) = hash_walk_through ($copy_H_category, [], \&deletekey);

                   if ($ref_array_v) {

                   print "$_ " foreach @$ref_key_list;
                   print " __ ";


                   my $total_cond = scalar(@$ref_array_v);

                    unless ($total_cond == 0) {


                         foreach my $uniq_interact (keys %{ $ref_h_all_corresp->{$Type_interaction} }) {

                           my $Count_Identical_Strictly      = 0;
                           my $Count_TypeConserved_Strictly  = 0;
                           my $Count_Identical_InTheArea     = 0;
                           my $Count_TypeConserved_InTheArea = 0;

                           foreach my $corresp_tocheck (@{ $ref_h_all_corresp->{$Type_interaction}{$uniq_interact} }) {
                                 #Identical_Strictly H1N1 1ru7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

                                 my $all_counts = conservation_amongst_conditions ($ref_array_v, $corresp_tocheck);

                                 my @counted = split (" ", $all_counts);
                                 $Count_Identical_Strictly += $counted[0];
                                 $Count_TypeConserved_Strictly += $counted[1];
                                 $Count_Identical_InTheArea += $counted[2];
                                 $Count_TypeConserved_InTheArea += $counted[3];

                                 }


                         my $Count_strictly_total = $Count_Identical_Strictly + $Count_TypeConserved_Strictly;

                         my $fairly_conserved   = 0.5 * $total_cond;     #At least 50% of concerned PDB (Host human for instance) are either Identical_Strictly or TypeConserved_Strictly
                         my $decently_conserved = 0.65 * $total_cond;
                         my $highly_conserved   = 0.8  * $total_cond;
                         my $mostly_conserved   = 0.9 * $total_cond;
                         my $fully_conserved    = 0.95 * $total_cond;
                         my $ratio_I_S   = $Count_Identical_Strictly/$total_cond;   my $ratio_TC_S   = $Count_TypeConserved_Strictly/$total_cond;
                         my $ratio_I_ITA = $Count_Identical_InTheArea/$total_cond;  my $ratio_TC_ITA = $Count_TypeConserved_InTheArea/$total_cond;
                         my $values_ratio = "$ratio_I_S $ratio_TC_S $ratio_I_ITA $ratio_TC_ITA";

                         if ($fairly_conserved <= $Count_strictly_total && $decently_conserved > $Count_strictly_total) {
                           $H_conservation{$flu}{$State}{$Type_interaction}{$uniq_interact}{"fairly_conserved"}{$ref_key_list->[0]}{$ref_key_list->[1]} = $values_ratio;
                           }
                             elsif ($decently_conserved <= $Count_strictly_total && $highly_conserved > $Count_strictly_total) {
                               $H_conservation{$flu}{$State}{$Type_interaction}{$uniq_interact}{"decently_conserved"}{$ref_key_list->[0]}{$ref_key_list->[1]} = $values_ratio;
                             }
                               elsif ($highly_conserved <= $Count_strictly_total && $mostly_conserved > $Count_strictly_total) {
                                   $H_conservation{$flu}{$State}{$Type_interaction}{$uniq_interact}{"highly_conserved"}{$ref_key_list->[0]}{$ref_key_list->[1]} = $values_ratio;
                                   }
                                   elsif ($mostly_conserved <= $Count_strictly_total && $fully_conserved > $Count_strictly_total) {
                                       $H_conservation{$flu}{$State}{$Type_interaction}{$uniq_interact}{"mostly_conserved"}{$ref_key_list->[0]}{$ref_key_list->[1]} = $values_ratio;
                                       }
                                       elsif ($fully_conserved  <= $Count_strictly_total) {
                                           $H_conservation{$flu}{$State}{$Type_interaction}{$uniq_interact}{"fully_conserved"}{$ref_key_list->[0]}{$ref_key_list->[1]} = $values_ratio;
                                           }
                       }

                     }
                   }
               }

                           print "\n--> DONE\n\n";

                           print "\n\nAssessing GLOBAL conservation  :\n";

                           assess_global_conservation (\@List_all_codePDB, $ref_h_all_corresp, \%H_conservation, $flu, $State, $Type_interaction, \@all_unique_conserved);


             undef $ref_h_all_corresp;
             }

                           print "--> DONE\n\n";

                           my $OutPath = "$path_environment/Interactions/Result/Conservation/$flu/$State/$Type_interaction";
                           h_json_and_dumper (\%H_conservation, $OutPath, "Conservation_values_$Type_interaction");
                           undef %H_conservation;
                           print uc ("$Type_interaction CONSERVATION VALUES DETERMINED FOR EVERY CLASSIFICATION\n\n");
                           print "...............................................................................................................................\n\n\n";
                           print "\n\n" x 5;


        }
      }
     }
    }
  }


           print "\n\n\n\n\n\n____ INDIVIDUAL INTERACTION CONSERVATION HAS BEEN ASSESSED  ____\n";
           print "____ INDIVIDUAL INTERACTION CONSERVATION HAS BEEN ASSESSED  ____\n";
           print "____ INDIVIDUAL INTERACTION CONSERVATION HAS BEEN ASSESSED  ____\n\n\n";
############## CONSERVATION PART END ###############
############## CONSERVATION PART END ###############
############## CONSERVATION PART END ###############



















  ####################   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   ####################
  ####################   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   ####################
  ####################   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   ####################
  ####################   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   ####################


  sub opendir_readdir {

        my $dir = shift;

        opendir(DIR, "$dir") or die "Unable to enter dir $dir:$!\n";
        my @FolderNames = readdir(DIR) or die "Unable to read $dir:$!\n";
        @FolderNames =  grep { ! m/^\./ } @FolderNames;
        close DIR;

        return @FolderNames;

  }






  sub get_all_corresponding_interaction_and_assess_similarity {
    my ($flu, $State, $Type_interaction, $Unique_interaction, $ref_H_state_subtype_code) = @_;

    my %temp_H;
    my @Split_uniq_interaction = split(" ", $Unique_interaction);

    foreach my $Subtype (keys %{ $ref_H_state_subtype_code->{$flu}{$State} }) {

      foreach my $CodePDB ( @{ $ref_H_state_subtype_code->{$flu}{$State}{$Subtype} } ) {


        my $Path = "$path_environment/Interactions/Result/$flu/$State/$Subtype/$CodePDB/";
        my $File = $Type_interaction . "_" . "$CodePDB.txt";
        open (FILE   , "<$Path$File") or die "Pb d'ouverture $Path$File : $!";


        while (<FILE>) {
          # H1N1 1RU7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299
          #  0    1      2        3   4   5  6    7       8        9     10   11 12  13  14     15     16      17

          chomp;
          my @Split = split(" ", $_);
          foreach (@Split) { $_ =~ s/^\s+|\s+$//g; }

          #Check whether HA1/HA2 location of interaction patners are the same between the two interactions
          #Check whether interaction type is the same between the two interactions
          my $HA1_HA2_OK_or_NOT = check_type_and_HA1_HA2 ($Split[6], $Split[13], $Split[2], \@Split_uniq_interaction);


          next if $HA1_HA2_OK_or_NOT ne "Ok";


          #Check whether residues are identical or not between the two interactions
          my $Conservation_residu = check_residu_conservation ($Split[4], $Split[11], \@Split_uniq_interaction);

          #Check whether numeration of each residues is the same between the two interactions
          my $Conservation_numerotation = check_numerotation_conservation ($Split[5], $Split[12], \@Split_uniq_interaction);



          #WE KEEP IF IDENTICAL residues OR mutation which allow to maintain the interaction type
          # AND  IF the numeration are Strictly the same or InTheArea (InTheArea usable for future work when we have more structure datas)
          # InTheArea conservation was not considered in the paper
          if ( $Conservation_residu       =~ m/Identical|TypeConserved/
          && $Conservation_numerotation =~ m/Strictly|InTheArea/) {

            my $Info = join ("_", $Conservation_residu, $Conservation_numerotation);
            my $ToPush_interaction = $Info . " " . $_;
            push @{ $temp_H{$Type_interaction}{$Unique_interaction} }, $ToPush_interaction;

          }
        }


        close FILE;
      }
     }
    return(\%temp_H);
  }





  sub check_type_and_HA1_HA2 {

       my ( $first, $second, $Type, $ref_SplitUniq) = @_;

       #Uniqinteraction :   $Splt[2], $Splt[4], $Splt[5], $Splt[6], $Splt[11], $Splt[12], $Splt[13]);
       #                     Type     Residue     Number     HA1      Residue    Number       HA1

        my $Info;

        if ($Type eq $ref_SplitUniq->[0]) {

            if ( ($ref_SplitUniq->[3] eq $first && $ref_SplitUniq->[6] eq $second)
              || ($ref_SplitUniq->[3] eq $second && $ref_SplitUniq->[6] eq $first) ) {

                  $Info = "Ok";
                }

                else { $Info = "NotOk"; }

              }

     return $Info;
   }





  sub check_residu_conservation {

     my ($first, $second, $ref_SplitUniq) = @_;
     #Uniqinteraction :   $Splt[2], $Splt[4], $Splt[5], $Splt[6], $Splt[11], $Splt[12], $Splt[13]);
     #                     Type     Residue     Number     HA1      Residue    Number       HA1

     my $Info;

     if ( ($ref_SplitUniq->[1] eq $first && $ref_SplitUniq->[4] eq $second)
       || ($ref_SplitUniq->[1] eq $second && $ref_SplitUniq->[4] eq $first) ) {

           $Info = "Identical";
         }

         else { $Info = "TypeConserved"; }

         return $Info;

  }




  sub check_numerotation_conservation {

     my ($first, $second, $ref_SplitUniq) = @_;

     my $Info;

     #VARIABLE FOR "InTheArea" Check. (InTheArea usable for future work when we have more structure datas)
     #CUTOFF difference of 5 positions max
     my $uniq_one_plus_five = $ref_SplitUniq->[2] + 5;    my $uniq_one_less_five = $ref_SplitUniq->[2] - 5;
     my $uniq_two_plus_five = $ref_SplitUniq->[5] + 5;    my $uniq_two_less_five = $ref_SplitUniq->[5] - 5;


     if (  ($first == $ref_SplitUniq->[2] && $second == $ref_SplitUniq->[5])
        || ($first == $ref_SplitUniq->[5] && $second == $ref_SplitUniq->[2]) ) {

          $Info = "Strictly";
        }


     elsif (  ( ($first >= $uniq_one_less_five && $first <= $uniq_one_plus_five) && ($second >= $uniq_two_less_five && $second <= $uniq_two_plus_five) )
           || ( ($first >= $uniq_two_less_five && $first <= $uniq_two_plus_five) && ($second >= $uniq_one_less_five && $second <= $uniq_one_plus_five) ) ) {

             $Info = "InTheArea";
           }

      else { $Info = "NotConserved"; }

         return $Info;

  }






  sub keep_mostsimilar_interaction {

  my ($ref_1, $ref_2, $ref_3, $ref_4, $interact) = @_;
  #my @Grep_1 = grep{ /Identical_Strictly/ }  @Grep;        my @Grep_2 = grep{ /TypeConserved_Strictly/ }  @Grep;
  #my @Grep_3 = grep{ /Identical_InTheArea/ } @Grep;        my @Grep_4 = grep{ /TypeConserved_InTheArea/ } @Grep;


  # Identical_Strictly H1N1 1RU7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299
  #        0            1    2        3       4  5   6   7    8       9      10     11   12 13 14    15     16    17       18


  my $Keep;

  my @first_tocheck  = ($ref_1, $ref_2);
  my @second_tocheck = ($ref_3, $ref_4);

  for my $ref (@first_tocheck) {

    if (scalar(@$ref) > 0) {

      my @temp_array = ();

        foreach my $interaction (@$ref) {

          my @Splt = split " ", $interaction; foreach (@Splt) { $_ =~ s/^\s+|\s+$//g; }
          # Identical_InTheArea H1N1 1RU7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

            my $first_res  = join(" ", $Splt[5], $Splt[6], $Splt[7]);
            my $second_res = join(" ", $Splt[12], $Splt[13], $Splt[14]);

            my $potential_unique_order   = "$first_res.+$second_res";
            my $potential_unique_reverse = "$second_res.+$first_res";

            unless ( grep{ $_ =~ /$potential_unique_order|$potential_unique_reverse/ } @temp_array ) {
                    push @temp_array, $interaction;
                  }
        }

        if (scalar(@temp_array) == 1 ) { $Keep = $temp_array[0]; return ($Keep, $ref); }
        else { my $ok; my $distance;
               foreach (@temp_array) {
                  my @split = split (" ", $_); foreach (@split) { $_ =~ s/^\s+|\s+$//g; }
                  if (!defined($distance)) { $ok = $_; $distance = $split[18]; }
                  elsif ($split[18] < $distance) { $ok = $_; $distance = $split[18]; }
                }
              return ($ok, $ref);
        }

   }
  }





    my @R_split = split (" ", $interact);
    #      Intra-Ionic  LYS  308 HA1  GLN 62 HA2
    foreach (@R_split) { $_ =~ s/^\s+|\s+$//g; }

    for my $ref (@second_tocheck) {
      #InTheArea so it could be more than 2 interactions *3 protomer retrieved : need to calculate to closest one
      if (scalar(@$ref) > 0) {
        my @temp_array = ();
        my $abs_value_sum;
        my @temp_array2 = ();

          foreach my $interaction (@$ref) {



            my @Splt = split " ", $interaction; foreach (@Splt) { $_ =~ s/^\s+|\s+$//g; }
            # Identical_InTheArea H1N1 1RU7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

            my $first_res  = join(" ", $Splt[5], $Splt[6], $Splt[7]);
            my $second_res = join(" ", $Splt[12], $Splt[13], $Splt[14]);

            my $potential_unique_order   = "$first_res.+$second_res";
            my $potential_unique_reverse = "$second_res.+$first_res";

                unless ( grep{ $_ =~ /$potential_unique_order|$potential_unique_reverse/ } @temp_array ) {

                    #Need the correspondance to calculate absolute values
                    my ($corresp_first, $corresp_second) = get_residu_inthearea_correspondance_to_unic_interaction ($R_split[2], $R_split[5], $Splt[6], $Splt[13]);

                    my $To_push = $corresp_first . "_" . $corresp_second . "__" . $interaction;
                    push @temp_array, $To_push;
                  }
                }


                foreach my $push (@temp_array) {
                  $push =~ m/^(\d+)_(\d+)__(.+)/;

                  my $abs_1 = abs($R_split[2] - $1);  my $abs_2 = abs($R_split[5] - $2); my $abssum = $abs_1 + $abs_2;
                  if ( ! defined($abs_value_sum) )   { $abs_value_sum = $abssum; push @temp_array2, $3; }
                  elsif ( $abssum <= $abs_value_sum ) { $abs_value_sum = $abssum; push @temp_array2, $3; }
                }

                if (scalar(@temp_array2) == 1 ) { $Keep = $temp_array2[0]; return ($Keep, $ref); }
                else {
                  my $ok; my $distance;
                       foreach (@temp_array2) {
                          my @split = split (" ", $_); foreach (@split) { $_ =~ s/^\s+|\s+$//g; }
                          if (!defined($distance)) { $ok = $_; $distance = $split[18]; }
                          elsif ($split[18] < $distance) { $ok = $_; $distance = $split[18]; }
                        }
                      return ($ok, $ref);

              }
            }

     }
}



  sub get_residu_inthearea_correspondance_to_unic_interaction {

  my ($Pos1, $Pos2, $Check1, $Check2) = @_;

  my $P1_plus5 = $Pos1 + 5; my $P1_less5 = $Pos1 - 5;
  my $P2_plus5 = $Pos2 + 5; my $P2_less5 = $Pos2 - 5;

  my $Ope_1 = $Pos1 - $Check1;
  my $Ope_2 = $Pos1 - $Check2;

  my $Ope_3 = $Pos2 - $Check1;
  my $Ope_4 = $Pos2 - $Check2;

  ### One of the 2 interacting residue is inthearea
  if ( ($Check1 >= $P1_less5 && $Check1 <= $P1_plus5) && !($Check2 >= $P1_less5 && $Check2 <= $P1_plus5) ) { return ($Check1, $Check2); }

  if ( ($Check2 >= $P1_less5 && $Check2 <= $P1_plus5) && !($Check1 >= $P1_less5 && $Check1 <= $P1_plus5) ) { return ($Check2, $Check1); }



  if ( ($Check1 >= $P2_less5 && $Check1 <= $P2_plus5) && !($Check2 >= $P2_less5 && $Check2 <= $P2_plus5) ) { return ($Check2, $Check1); }

  if ( ($Check2 >= $P2_less5 && $Check2 <= $P2_plus5) && !($Check1 >= $P2_less5 && $Check1 <= $P2_plus5) ) { return ($Check1, $Check2); }



  ### The 2 interacting residue are inthearea
  if (  ($Check1 >= $P1_less5 && $Check1 <= $P1_plus5) && ($Check2 >= $P1_less5 && $Check2 <= $P1_plus5)  ) {
            if ( abs($Ope_3) >= abs($Ope_4) )    { return ($Check1, $Check2); }
            elsif ( abs($Ope_3) <= abs($Ope_4) ) { return ($Check2, $Check1); }
  }



  if (  ($Check1 >= $P2_less5 && $Check1 <= $P2_plus5) && ($Check2 >= $P2_less5 && $Check2 <= $P2_plus5)  ) {
            if ( abs($Ope_1) >= abs($Ope_2) )    { return ($Check2, $Check1); }
            elsif ( abs($Ope_1) <= abs($Ope_2) ) { return ($Check1, $Check2); }
  }

  return;
  }







  sub interact_symmetry_check {

  my ($Keep, $ref_grep_keep, $ref_grep_base) = @_;
  # $Keep
  #TypeConserved_InTheArea H18N11	4mc5	Inter-Hydrophobic	CE2 TYR 458 HA2 Proto1 old=423   |     CD1 LEU 463 HA2 Proto3 old=428 ->   3.96914474414828
  #             0            1      2            3          4  5   6   7   8      9        10    11   12  13  14   15      16    17       18


  my @Split = split (" ", $Keep);
  foreach (@Split) { $_ =~ s/^\s+|\s+$//g; }

      my $regex_keep = "$Split[5]\\s$Split[6]\\s$Split[7].+$Split[12]\\s$Split[13]\\s$Split[14]";
      my $regex_keep_reverse = "$Split[12]\\s$Split[13]\\s$Split[14].+$Split[5]\\s$Split[6]\\s$Split[7]";

      my @filter_interact = grep{ /$regex_keep|$regex_keep_reverse/ }  @$ref_grep_keep;
      #Grep all interactions corresponding to $Keep


      my $regex_proto_nb1;
      my $regex_proto_nb2;
      my $regex_proto_nb3;

      if ($Split[3] =~ /Intra/) {
        $regex_proto_nb1 = "Proto1";
        $regex_proto_nb2 = "Proto2";
        $regex_proto_nb3 = "Proto3";
       }

      elsif ($Split[3] =~ /Inter/) {
        $regex_proto_nb1 = "Proto1.+Proto2|Proto2.+Proto1";
        $regex_proto_nb2 = "Proto1.+Proto3|Proto3.+Proto1";
        $regex_proto_nb3 = "Proto2.+Proto3|Proto3.+Proto2";
       }


            my @grep_5 = grep{ /$regex_proto_nb1/ }  @filter_interact;
            my @grep_6 = grep{ /$regex_proto_nb2/ }  @filter_interact;
            my @grep_7 = grep{ /$regex_proto_nb3/ }  @filter_interact;
            #Grep l'interaction pour chancun des trois proto qu'on soit en intra ou en interproto


            if (@grep_5 && @grep_6 && @grep_7) {
              return $Keep;
             }

            else {
              my @array_ref_grep567 = (\@grep_5, \@grep_6, \@grep_7);
              my @array_grep_to_splice = ($ref_grep_keep, $ref_grep_base);

              foreach my $refgrep567 (@array_ref_grep567) {
              if (@$refgrep567) {
                foreach my $interact (@$refgrep567) {


                    foreach my $grep_splice (@array_grep_to_splice) {

                      my $i = 0;
                      while ($i <= $#$grep_splice) {

                         splice (@$grep_splice, $i, 1) if $grep_splice->[$i] eq $interact;
                         $i++;
                       }

                    }

                 }
                }
               }
            return;
            }

    }











  sub splice_from_array {

  my ($Array_to_splice_in, $Motif_to_check, $CodePDB) = @_;

  if (ref($Motif_to_check) eq "ARRAY") {
                            foreach my $motif (@$Motif_to_check ) {
                            my $i = 0;
                                foreach (@$Array_to_splice_in) {  $i++;
                                                                  my @splice = splice(@$Array_to_splice_in, $i, 1) if $_ eq $motif;
                                                                }
                          }
    }


   elsif (ref($Motif_to_check) eq "SCALAR") {
                              my $i = 0;
                              foreach (@$Array_to_splice_in) {
                                  $i++;
                                                                my @splice = splice(@$Array_to_splice_in, $i, 1) if $_ eq $Motif_to_check;
                                                                }
    }
  return;
    }





    sub h_json_and_dumper {
          my ($H_RealInteractions, $Outpath, $Info) = @_;

          my $path_json = "$Outpath/JSON"; my $path_dump = "$Outpath/DUMPER";
          make_path($path_json);           make_path($path_dump);

          open(HFILE,">$path_json/$Info.dat") or die "Impossible to open : $path_json/$Info.dat : $!";;
          my $path_json_file = "$path_json/$Info.dat";
          my $JSON_H = encode_json($H_RealInteractions);
          print HFILE $JSON_H;
          close HFILE;

          print "\n------> Hash $Info : Saved\n";

          open(HASHDUMP,">$path_dump/$Info.txt");
          print HASHDUMP Dumper ($H_RealInteractions);
          close HASHDUMP;

      return $path_json_file;
                             }



      sub get_json_data_hardsaved {

      my $path = shift;

      open(DATFILE,"<$path") or warn "Impossible to open : $path --> $!";

      my $JSONdata = <DATFILE>;

      my $h_ref = decode_json($JSONdata);

      close DATFILE;

      return $h_ref;
      }





      sub hash_walk_through {
          my ($hash, $key_list, $callback) = @_;

          my @keys = qw(SubtypesSameHA Subtypes Hosts HA_Clades HA_Groups Globally);
          my $KEY; my $VALUE; my $ref_KEY_LIST;

          while (my ($key, $value) = each %$hash) {
              # Keep track of the hierarchy of keys, just in case
              # our callback needs it.
              push @$key_list, $key;
              if ($key_list->[0] ne $key) {
              if ( (any {$key_list->[0] eq $_} (@keys)) && (any {$key eq $_} (@keys)) ) { shift @$key_list;  }
                  }
      #        print "xxxxxxxxxx $key $value  \n";
      #        print Dumper $key_list;
      #      print "\n";
      #print Dumper $key_list;
              if (ref($value) eq 'HASH') {

              if (!%$value) { delete $hash->{$key}; next;}
                  # Recurse.
                  #GOING INTO THE SECOND CALL BY FINISHING THE FIRST ONE : VERY IMPORTANT FOR RECURSION
                  return hash_walk_through($value, $key_list, $callback);
              }

              elsif (defined($KEY) && defined($VALUE)) { return ($KEY, $VALUE, $ref_KEY_LIST); }

              elsif (ref($value) eq 'ARRAY') {
                  # Otherwise, invoke our callback, passing it
               # the current key and value, along with the
                  # full parentage of that key.
                 $KEY = $key; $VALUE = $value; $ref_KEY_LIST = $key_list;

                  #Callback other function
                  $callback->($key, $hash);

                  return($KEY, $VALUE, $ref_KEY_LIST);
              }
              if (defined($KEY) && defined($VALUE)) { return ($KEY, $VALUE, $ref_KEY_LIST); }
        }
      }


      sub deletekey {
      my ($key, $ref_hash) = @_;


      if (exists($ref_hash->{$key})) { delete $ref_hash->{$key};}

      else {  print "no key ....";}

      return;

      }




      sub conservation_amongst_conditions {
        my ($ref_array_v, $corresp_tocheck) = @_;

        my $Count_Identical_Strictly = 0;
        my $Count_TypeConserved_Strictly = 0;
        my $Count_Identical_InTheArea = 0;
        my $Count_TypeConserved_InTheArea = 0;

                  my @Split = split (" ", $corresp_tocheck);
                  foreach (@Split) { $_ =~ s/^\s+|\s+$//g; }
                  #Identical_Strictly H1N1 1ru7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

                  my $ToCheck = uc ($Split[2]);

                  my @grep = grep { /$ToCheck|$Split[2]/ } @$ref_array_v;


                  if (@grep != 0) {
                   $Count_Identical_Strictly++       if $corresp_tocheck =~ /Identical_Strictly/;
                   $Count_TypeConserved_Strictly++   if $corresp_tocheck =~ /TypeConserved_Strictly/;
                   $Count_Identical_InTheArea++      if $corresp_tocheck =~ /Identical_InTheArea/;
                   $Count_TypeConserved_InTheArea++  if $corresp_tocheck =~ /TypeConserved_InTheArea/;
                        }

            my $counts = join(" ", $Count_Identical_Strictly, $Count_TypeConserved_Strictly, $Count_Identical_InTheArea, $Count_TypeConserved_InTheArea);

            return($counts);
            }






            sub assess_global_conservation {
              my ($ref_array_v, $ref_all_corresponding_interactions, $ref_H_conservation, $flu, $state, $type, $all_unique_conserved) = @_;

              my $total_cond = scalar(@$ref_array_v);
              my $fairly_conserved   = 0.5 * $total_cond;     #At least 50% of concerned PDB (Host human for instance) are either Identical_Strictly or TypeConserved_Strictly
              my $decently_conserved = 0.65 * $total_cond;
              my $highly_conserved   = 0.8  * $total_cond;
              my $mostly_conserved   = 0.9 * $total_cond;
              my $fully_conserved    = 0.95 * $total_cond;


               for my $uniq_interact (keys (%{ $ref_all_corresponding_interactions->{$type} })) {
                my $Count_Identical_Strictly      = 0;
                my $Count_TypeConserved_Strictly  = 0;
                my $Count_Identical_InTheArea     = 0;
                my $Count_TypeConserved_InTheArea = 0;




                  foreach my $Interact ( @{ $ref_all_corresponding_interactions->{$type}{$uniq_interact} } ) {
                        #Identical_Strictly H1N1 1ru7 Intra-Ionic NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

                        my @Split = split (" ", $Interact);
                        foreach (@Split) { $Interact =~ s/^\s+|\s+$//g; }

                        my $ToCheck = uc($Split[2]);
                        my @grep = grep { /$ToCheck|$Split[2]/ } @$ref_array_v;


                        if (@grep != 0) {
                         $Count_Identical_Strictly++      if $Interact =~ /Identical_Strictly/;
                         $Count_TypeConserved_Strictly++  if $Interact =~ /TypeConserved_Strictly/;
                         $Count_Identical_InTheArea++     if $Interact =~ /Identical_InTheArea/;
                         $Count_TypeConserved_InTheArea++ if $Interact =~ /TypeConserved_InTheArea/;
                         }
                    }


                    my $Count_strictly_total = $Count_Identical_Strictly + $Count_TypeConserved_Strictly;
                    my $ratio_I_S   = $Count_Identical_Strictly/$total_cond;   my $ratio_TC_S   = $Count_TypeConserved_Strictly/$total_cond;
                    my $ratio_I_ITA = $Count_Identical_InTheArea/$total_cond;  my $ratio_TC_ITA = $Count_TypeConserved_InTheArea/$total_cond;
                    my $values_ratio = "$ratio_I_S $ratio_TC_S $ratio_I_ITA $ratio_TC_ITA";

                    if ($fairly_conserved <= $Count_strictly_total && $decently_conserved > $Count_strictly_total) {
                    $ref_H_conservation->{$flu}{$state}{$type}{$uniq_interact}{"fairly_conserved"}{"Globally"}{"Globally"} = $values_ratio;
                    }
                      elsif ($decently_conserved <= $Count_strictly_total && $highly_conserved > $Count_strictly_total) {
                        $ref_H_conservation->{$flu}{$state}{$type}{$uniq_interact}{"decently_conserved"}{"Globally"}{"Globally"} = $values_ratio;
                      }
                        elsif ($highly_conserved <= $Count_strictly_total && $mostly_conserved > $Count_strictly_total) {
                            $ref_H_conservation->{$flu}{$state}{$type}{$uniq_interact}{"highly_conserved"}{"Globally"}{"Globally"} = $values_ratio;
                            }
                            elsif ($mostly_conserved <= $Count_strictly_total && $fully_conserved > $Count_strictly_total) {
                                $ref_H_conservation->{$flu}{$state}{$type}{$uniq_interact}{"mostly_conserved"}{"Globally"}{"Globally"} = $values_ratio;
                                }
                                elsif ($fully_conserved  <= $Count_strictly_total) {
                                    $ref_H_conservation->{$flu}{$state}{$type}{$uniq_interact}{"fully_conserved"}{"Globally"}{"Globally"} = $values_ratio;
                                    }
                      push @$all_unique_conserved, $uniq_interact if $Count_strictly_total > 0.6;
                                # H_conservation : key1 = $flu, key2 = $state, key3 = (conservation : 75, 85 or 95), key4 = uniq interaction, $key5 = Folder_group, $key6 = % conservation
                  }

            return;
            }







            sub get_groups {

            my ($interaction, $ref_all_uniq, $ref_all_groups) = @_;

            my $Count = 0;

            $interaction =~ m/.+-.+\s(.+\s\d+\sHA\d)\s(.+\s\d+\sHA\d)/;
            #$Splt[2], $Splt[4], $Splt[5], $Splt[6], $Splt[11], $Splt[12], $Splt[13]
            # Intra-Ionic LYS 308 HA1   GLN 62 HA2

            my @Split_1 = split (" ", $interaction);
            my @Split_2 = split ("-", $Split_1[0]);
            return if $Split_2[1] eq "Hydrogen";
            return if $Split_2[1] eq "Repulsive";

            my @first = grep{ /$1/ } @$ref_all_uniq;
            my @second = grep{ /$2/ } @$ref_all_uniq;

            my @group;
            push @group, $interaction;

            my $ref_first = \@first;
            my $ref_second = \@second;

            foreach my $ref_grep ($ref_first, $ref_second) {
                LOOPY:
                foreach (@$ref_grep) {
                  my @Split_1 = split (" ", $_);
                  my @Split_2 = split ("-", $Split_1[0]);
                  next LOOPY if $Split_2[1] eq "Hydrogen";
                  next LOOPY if $Split_2[1] eq "Repulsive";
                  push @group, $_;
                }
              }

            @group = uniq @group;
            return if scalar(@group) <= 1;

            push @$ref_all_groups, \@group;

            return;
            }




            sub keep_best_interact_if_same_but_different_conservation {
            #Meaning of this sub is to check if there is already a similar interaction and keep the best one
            #key H18N11 : Intra-Hydrogen_SD_C GLY 41 HA1 THR 351 HA1' => '0 1 0 0',
            #key H18N11 : Intra-Hydrogen_SD_C SER 41 HA1 THR 351 HA1' => '1 0 0 0'

            my ($interact, $ref_value, $ref_temp_H) = @_;


            my @split = split (" ", $interact);
            my @split_v = split (" ", $ref_value);

            my @grep = grep { $_ =~ m/$split[0].+$split[2].+$split[5]|$split[0].+$split[5].+$split[2]/ } (keys %$ref_temp_H);

            if (@grep) {
              my $to_keep;
              my $value_to_keep;

                for (@grep) {
                  my @splt_val = split (" ", $ref_temp_H->{$_});
                  $to_keep = $_ if $splt_val[0] > $split_v[0];
                  $to_keep = $_ if $splt_val[0] == $split_v[0] && $splt_val[1] > $split_v[1];
                  $to_keep = $_ if $splt_val[0] == $split_v[0] && $splt_val[1] == $split_v[1] && $splt_val[2] > $split_v[2];
                  $to_keep = $_ if $splt_val[0] == $split_v[0] && $splt_val[1] == $split_v[1] && $splt_val[2] == $split_v[2] && $splt_val[3] > $split_v[3];

                  }

              #  $value_to_keep = $ref_temp_H->{$to_keep} if defined($to_keep);
              return if defined($to_keep);

              if (!defined($to_keep)) {
                                        foreach my $key (keys %$ref_temp_H) {
                                            delete $ref_temp_H->{$key} if any{$key eq $_} (@grep);
                                          }

                                          $value_to_keep = $ref_value;
                                          $to_keep = $interact;

                                          $ref_temp_H->{$interact} = $ref_value;
                                          return;
                                      }

              }

              else {
                $ref_temp_H->{$interact} = $ref_value;
                return;
                }

            }
