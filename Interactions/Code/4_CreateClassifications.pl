use warnings;
use strict;

use File::Path qw( make_path );
use List::Util qw(sum uniq any max reduce);
use JSON;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;




    #Check DUPLICATED PDB CODE AND CREATED GROUPS OF STRUCTURES (WITH ALL ASSOCIATED INFORMATIONS such as structure resolution) DEPENDING ON HA CLASSIFICATION : Group/clade/subtypeSameHA/Subtype/Host...


    my $path_environment = "/Users/valentin.ozl/GitHub";
    #    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment


###########################
#CHECK IF WE CLEANED ALL STRUCTURES WHICH ARE REPRESENTED IN OUR DATASET (Extended Fig 1)

    open (TABLEAUPDB, "<$path_environment/Interactions/Data/TableauPDB.csv") or die "Unable to open file TableauPDB.csv :$!\n";
    #Data set with all structural informations (Extended Data Fig 1)

    my @ListCodePDBinTable = ();

    my $CountLine = 0;

    my @Subtypes = ();
    my @Hosts = ();
    my @TableauPDB = ();
    my @TableauPDB_withRefPaper = ();
    my $Count = 0;



    while ( my $Line = <TABLEAUPDB>) {

      $CountLine++;

      next if $CountLine < 10; #Headers

      if ($Line =~ /^(([[:alnum:]]){4})/) {

          push @ListCodePDBinTable, $1;

          #EXEMPLE :
          #3UBQ;H1N1;A/California/04/2009;Human;G205C/R220C (HA1);3SLN;X-RAY DIFFRACTION;2;Structural Characterization of the Hemagglutinin Receptor Specificity from the 2009 H1N1 Influenza Pandemic.;;;;;;;;;;;;;;;;
          my ($PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution, $Paper) = split (";", $Line);

          $Line = join ("|", $PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution);
          my $LineWithPaper = join ("|", $PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution, $Paper);

          #Create lists with all PDB CODE informations
          push @TableauPDB, $Line; push @TableauPDB_withRefPaper, $LineWithPaper;
          push @Subtypes, $Subtype;
          push @Hosts, $Host;

        }
    }

    #Put all infos about PDB code in JSON FILE
    open(ARRAYCLEAN,">", "$path_environment/Interactions/Data/ArrayCleanWithRefTableauPDB.dat");
    my $JSON = encode_json(\@TableauPDB_withRefPaper);   print ARRAYCLEAN $JSON;




    my $Warning =0;
    my @UniqListCodePDBinTable = uniq @ListCodePDBinTable;
      foreach my $Uniq (@UniqListCodePDBinTable) {

        my @Grep = grep (m/^$Uniq/, @ListCodePDBinTable);
        my $Number = scalar (@Grep);
        $Warning++, print "warning doublon $Uniq\n" if scalar(@Grep) >= 2;
       }

      print "\nTableauPDB.csv OK --> we can procede FastaFile\n" if $Warning == 0;
      close TABLEAUPDB;
      close ARRAYCLEAN;

############################









#CLASS STRUCTURES ACCORDING TO CLASSIFICATIONS (Exemple: All H1N1 structures, all Clade 3 structures...)


# PREEXISTING FILES CLASSED BY STATE OF HA : Prefusion / Intermediary/Postfus / FluB
my %FilesToProcede = ( "Prefusion" => "$path_environment/Interactions/Data/ArrayTableauPDB/FluA/Prefusion",
                       "Postfusion_Intermediary" => "$path_environment/Interactions/Data/ArrayTableauPDB/FluA/Postfusion_Intermediary",
                       "PrefusionFluB" => "$path_environment/Interactions/Data/ArrayTableauPDB/FluB/Prefusion",
                       );


while ((my $ConditionHA, my $Directory) = each %FilesToProcede)
        {

          $ConditionHA = "Prefusion" if $ConditionHA =~ m/Prefusion/;

          #CSV files of structural infos (FluA prefusion, FluA postfusion/intermediates or FluB)
          open (TABLEPDB, "<", "$Directory/TableauPDB_$ConditionHA.csv") or die "Unable to open file TableauPDB_$ConditionHA.csv :$!\n";

          @ListCodePDBinTable = ();
          $CountLine = 0;
          @Subtypes = ();
          @Hosts = ();
          @TableauPDB = ();
          @TableauPDB_withRefPaper = ();



          while ( my $Line = <TABLEPDB>) {

            $CountLine++;
            #Headers
            next if $CountLine < 10;
            if ($Line =~ /^(([[:alnum:]]){4})/) {

                push @ListCodePDBinTable, $1;

                #3UBQ;H1N1;A/California/04/2009;Human;G205C/R220C (HA1);3SLN;X-RAY DIFFRACTION;2;Structural Characterization of the Hemagglutinin Receptor Specificity from the 2009 H1N1 Influenza Pandemic.;;;;;;;;;;;;;;;;
                my ($PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution, $Paper) = split (";", $Line);

                $Line = join ("|", $PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution);
                my $LineWithPaper = join ("|", $PDB, $Subtype, $Strain, $Host, $Mutation, $Ligand, $Method, $Resolution, $Paper);


                #Create lists with all PDB CODE informations FOR THIS SPECIFIC STATE
                push @TableauPDB, $Line;      push @TableauPDB_withRefPaper, $LineWithPaper;
                push @Subtypes, $Subtype;
                push @Hosts, $Host;

              }
          }




          my %HashCount = ();
          my %HashGroupsSameHA = ();
          @Subtypes = uniq @Subtypes;
          @Hosts = uniq @Hosts;

          #DETAILED CLASSIFICATIONS : SUBTYPES_sameHA/Subtypes/GROUP/CLADE/HOST
          #Subtypes same HA exemples : H1, H6, H7..
          #Subtypes exemples : H1N1, H6N8..
          my @Group1 = ("H1", "H2", "H5", "H6", "H8", "H9", "H11", "H12", "H13", "H16", "H17", "H18");
          my @Group2 = ("H3", "H4", "H7", "H10", "H14", "H15");

          my @CladeH1_1 = ("H1", "H2", "H5", "H6");
          my @CladeH9_1 = ("H8", "H9","H12");
          my @CladeH11_1 = ("H11","H13", "H16");
          my @CladeH17_1 = ("H17", "H18");

          my @CladeH3_2 = ("H3", "H4","H14");
          my @CladeH7_2 = ("H7", "H10", "H15");

          my %HashAllInfos = ();

          my %ToPassToSub = ( "Subtypes"     => \@Subtypes,

                              "HA_Groups"    => {
                                          "Group1" => \@Group1,
                                          "Group2" => \@Group2,
                                    },

                              "HA_Clades"    => {
                                          "CladeH1_1" => \@CladeH1_1,
                                          "CladeH9_1" => \@CladeH9_1,
                                          "CladeH11_1" => \@CladeH11_1,
                                          "CladeH17_1" => \@CladeH17_1,
                                          "CladeH3_2" => \@CladeH3_2,
                                          "CladeH7_2" => \@CladeH7_2,
                                    },

                              "Hosts"         => \@Hosts,

                                    );


              # CREATE SubtypeSameHA classification AND ADD IN Hash %ToPassSub
              for (my $i = 1; $i <= 18; $i++) {

                my $GroupOfSubtypesSameHA = "H$i";
                my @Group = grep { m/($GroupOfSubtypesSameHA)N/ } @Subtypes;
                $ToPassToSub{"SubtypesSameHA"}->{"H$i"} = \@Group;
              }

              # CREATE Subtype/host classification AND ADD IN Hash %ToPassSub
              foreach $_ (@Subtypes){
                my @Hostts = ();
                foreach my $Hostt (@Hosts){
                  push @Hostts, $Hostt;
                }
                $ToPassToSub{"Subtypes_Hosts"}->{$_} = \@Hostts;
              }







              #CREATE HASH WITH ALL CLASSIFICATION INFORMATION AND CORRESPONDING PDB CODES
              while ((my $KeyName1, my $Two) = each %ToPassToSub)
              {

                if (ref $Two eq "ARRAY") { Makehashallinfos (\@TableauPDB, \%HashAllInfos, \%HashCount, $KeyName1, undef, $Two); next;}


                elsif (ref $Two eq "HASH") {
                  while ((my $SecondKeyOrArray, my $Three) = each %$Two)
                  {
                    if (ref $SecondKeyOrArray eq "" && ref $Three eq "ARRAY") { Makehashallinfos (\@TableauPDB, \%HashAllInfos, \%HashCount, $KeyName1, $SecondKeyOrArray, $Three);}
                  }
                    }
                  }







                  # ALL STRUCTURES INFOS OF EXTENDED DATA FIG 1 AS AN HASH (IF NEEDEED)
                  open(ARRAYFILE,">", "$Directory/FileArrayTableauPDB_WithoutRefPaper_$ConditionHA.dat");
                  my $JSONdata1 = encode_json(\@TableauPDB);   print ARRAYFILE $JSONdata1;

                  open(ARRAYFILE2,">", "$Directory/ArrayTableauPDB_WithoutRefPaper_$ConditionHA.txt");
                  print ARRAYFILE2 Dumper ( \@TableauPDB );

                  open(ARRAYFILE3,">", "$Directory/FileArrayTableauPDB_WithRefPaper_$ConditionHA.dat");
                  my $JSONdata3 = encode_json(\@TableauPDB_withRefPaper);   print ARRAYFILE3 $JSONdata3;




                  # HASH OF DETAILED CLASSIFICATIONS (IF NEEDED)
                  open(HASHFILE1,">", "$Directory/HashAllConditions_$ConditionHA.dat");
                  my $JSONdata4 = encode_json(\%ToPassToSub);   print HASHFILE1 $JSONdata4;

                  open(HASHFILE1BIS,">", "$Directory/HashAllConditions_$ConditionHA.txt");
                  print HASHFILE1BIS Dumper (\%ToPassToSub);



                  # HASH OF CLASSIFIED STRUCTURES ACCORDING TO CLASSIFICATIONS
                  open(HASHFILE2,">", "$Directory/HashAllConditionsGrep_$ConditionHA.dat");
                  my $JSONdata5 = encode_json(\%HashAllInfos);   print HASHFILE2 $JSONdata5;

                  open(HASHFILE2BIS,">", "$Directory/HashAllConditionsGrep_$ConditionHA.txt");
                  print HASHFILE2BIS Dumper ( \%HashAllInfos );


                  # HASH OF NUMBER OF STRUCTURES PER CLASSIFICATION
                  open(HASHFILE3,">", "$Directory/HashAllConditionsCount_$ConditionHA.dat");
                  my $JSONdata6 = encode_json(\%HashCount);   print HASHFILE3 $JSONdata6;

                  open(HASHFILE3BIS,">", "$Directory/HashAllConditionsCount_$ConditionHA.txt");
                  print HASHFILE3BIS Dumper ( \%HashCount );



                  close TABLEPDB;
                  close ARRAYFILE;
                  close ARRAYFILE2;
                  close ARRAYFILE3;
                  close HASHFILE1;
                  close HASHFILE1BIS;
                  close HASHFILE2;
                  close HASHFILE2BIS;
                  close HASHFILE3;
                  close HASHFILE3BIS;

                  print "\nFiles ArrayTableauPDB done and saved : $ConditionHA condition\n\n"

  }




########################################################################################################
########## SUB ##########



             sub Makehashallinfos {
               my ($refTableauPDB, $refHashAllInfos, $refHashCount, $KeyName, $KeyName2, $refArray1) = @_;


                  #If non nested hash part
                  if (defined $KeyName && ! defined $KeyName2) {

                    foreach my $ToGrep ( @{ $refArray1 } ) {

                      #Grep subtypes to get all codePDB belonging to this subtype
                      my @Grep = grep { m/$ToGrep/ } @{ $refTableauPDB };
                      my $Count = scalar (@Grep);

                      # CREATE HASH WITH Classification info (Eg Subtypes), Precise Classification (Eg H1N1), and list containing all corresponding pdb codes + informations
                      $refHashAllInfos->{$KeyName}{$ToGrep} = [@Grep];
                      ### If we use \@Grep instead, it will change in the hash if we change @Grep. [@Grep] makes a ref with a copy of @Grep
                      $refHashCount->{$KeyName}{$ToGrep} = $Count;
                                                        }
                                                      }




                  #We don't have enough data for Subtype_Hosts at the moment
                  #If Nested hash part
                  elsif (defined $KeyName && defined $KeyName2 && $KeyName ne "Subtypes_Hosts") {

                     my @Grep = ();
                     my $Count = 0;

                              foreach my $ToGrep ( @{ $refArray1 } ) {

                                my @ProcedeGrep = ();

                                #Process for SubtypeSameHA, HA_Clades and HA_Groups
                                @ProcedeGrep = grep { m/($ToGrep)/ } @{ $refTableauPDB } if $KeyName eq "SubtypesSameHA";
                                @ProcedeGrep = grep { m/($ToGrep)N/ } @{ $refTableauPDB } if $KeyName eq "HA_Clades" || $KeyName eq "HA_Groups";
                                my $ProcedeCount = scalar (@ProcedeGrep);

                                push (@Grep, @ProcedeGrep);
                                $Count += $ProcedeCount;

                                                  }

                      # CREATE HASH WITH Classification info (Eg HA_Clades or HA_Groups for instance), Precise Classification (Eg CladeH1_1 or Group1), and list containing all corresponding pdb codes + informations
                      $refHashAllInfos->{$KeyName}{$KeyName2} = [@Grep];
                      ### If we use \@Grep instead, it will change in the hash if we change @Grep. [@Grep] makes a ref with a copy of @Grep
                      $refHashCount->{$KeyName}{$KeyName2} = $Count;
                                                  }







                        elsif (defined $KeyName && defined $KeyName2 && $KeyName eq "Subtypes_Hosts") {



                                  foreach my $ToGrep ( @{ $refArray1 } ) {
                                  my @Grep = grep { m/($KeyName2)/ && m/$ToGrep/ } @{ $refTableauPDB };
                                  my $Count = scalar (@Grep);


                                  unless (scalar (@Grep) == 0) {
                                  $refHashAllInfos->{$KeyName}{$KeyName2}{$ToGrep} = [@Grep];
                                  $refHashCount->{$KeyName}{$KeyName2}{$ToGrep} = $Count;
                                  }
                                                                                                  }
                                                                                  }
                            }
