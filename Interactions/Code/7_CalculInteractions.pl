use warnings;
use strict;

use File::Path qw(make_path);

use List::MoreUtils qw(any);

use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use JSON;

use Sort::Naturally;

use Parallel::ForkManager;
#███████████████████






########### COMPUTATION OF AA INTERACTIONS AT ATOMIC SCALE



my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



######################################################## INITIALISATION ############################################################################

my $info = Sys::Info->new;
my $cpu  = $info->device( 'CPU' );
my $NumberCPU = $cpu->count;
my $NumberCPUx2 = $NumberCPU *2;
my $forks = $NumberCPUx2 or die "Usage: $0 N\n";

my %Hash_Interactions_CorrespondingAtoms = (
"AtomIonicPositive" => ["NE", "NH1", "NH2", "NZ", "NE2", "ND1"],
"AtomIonicNegative" => ["OD1", "OD2", "OE1", "OE2"],
"AtomDonorProton" => ["OG", "OG1", "NE2", "ND2", "ND1", "NE2", "NZ", "NE", "NH1", "NH2", "OH", "NE1"],  #Sauf N backbone
"AtomAcceptorProton" => ["OG", "OG1", "OE1", "OE2", "OD1", "OD2", "ND1", "NE2", "OH"],  #Sauf O backbone
"AtomHydrophobic" => ["CB", "CG", "CE", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "CE1", "CZ", "CG1", "CG2", "CD", "CH2"]
);

my %Hash_Interactions_CorrespondingResidues = (
"ResiduIonicPositive" => ["ARG", "LYS", "HIS"],
"ResiduIonicNegative" => ["ASP", "GLU"],
"ResiduIonic"        => ["ARG", "LYS", "HIS", "ASP", "GLU"],
"ResiduDonorProton"  => ["SER", "THR", "GLN", "ASN", "HIS", "LYS", "ARG", "TYR", "TRP"],   #Sauf N donneur backbone de tous les AA
"ResiduAcceptorProton" => ["SER", "THR", "GLU", "ASP", "GLN", "ASN", "HIS", "TYR"],  #Sauf O accepteur backbone de tous les AA
"ResiduHydrophobic" => ["ALA", "MET", "TRP", "PHE", "TYR", "VAL", "LEU", "ILE", "PRO"]
);


my @AllInteractionTypes = ("All", "Ionic", "Repulsive", "RepulsivePositive", "RepulsiveNegative", "Hydrogen", "Hydrogen_SD_SD", "Hydrogen_SD_C", "Hydrophobic");


my @Indications = ( "Total", "Intra", "Inter", "HA1-HA1", "HA1-HA2", "HA2-HA2",
                      "IntraHA1-HA1", "IntraHA1-HA2", "IntraHA2-HA2",
                      "InterHA1-HA1", "InterHA1-HA2", "InterHA2-HA2" );

my $CountPDB = 0;

print "...\nWill process all following interaction types :\n";
print "$_ ...\n" foreach @AllInteractionTypes;
print "\n.........\nChange the array named AllInteractionsTypes if we add new interaction types\n.........\n\n";
print "\n.........\nMachine CPUs = $NumberCPU  /  Will use $forks fork for parallel code\n.........\n\n";


my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");


######################################################## INITIALISATION ############################################################################




my $time_beginning = gmtime();
print "\n\n\n\n__________ $time_beginning : Code start time\n\n\n\n\n\n";





foreach my $flu (@Flu) {


  foreach my $StructureState (@StructureStates){

    unless ($flu eq "FluB" and $StructureState eq "Intermediates" or $flu eq "FluB" and $StructureState eq "Postfusion"){

    my $dirInputSubtype = "$path_environment/Interactions/Data/PDB_Clean_aligned/Trimeric/$flu/$StructureState";
    my @SubtypesInputToProcess = opendir_readdir ($dirInputSubtype);



        foreach my $Subtype (@SubtypesInputToProcess) {

          next if ($Subtype =~ m/^\./);
          my $PathToCodesPDB = $dirInputSubtype . "/" . $Subtype;
          my @CodesPDBInputToProcess = opendir_readdir ($PathToCodesPDB);

            foreach $_ (@CodesPDBInputToProcess) {


              next if ($_ =~ m/^\./);
              $_ =~ m/PDB_Clean_aligned_(.+)\.txt/;   my $PathFileInput = $PathToCodesPDB . "/" . $_;
              my $CodePDB = $1;



              my $OutPath = "$path_environment/Interactions/Result/$flu/$StructureState/$Subtype/$CodePDB";
              #Create directories for the Output files
              make_path("$OutPath");

              #OPEN OUTPUT files
              open (INTERIONIC   , ">$OutPath/Inter-Ionic_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAIONIC  , ">$OutPath/Intra-Ionic_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERREPULSIVE   , ">$OutPath/Inter-Repulsive_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAREPULSIVE  , ">$OutPath/Intra-Repulsive_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERREPULSIVEPOSITIVE   , ">$OutPath/Inter-RepulsivePositive_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAREPULSIVEPOSITIVE  , ">$OutPath/Intra-RepulsivePositive_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERREPULSIVENEGATIVE   , ">$OutPath/Inter-RepulsiveNegative_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAREPULSIVENEGATIVE   , ">$OutPath/Intra-RepulsiveNegative_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERHYDROGEN   , ">$OutPath/Inter-Hydrogen_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAHYDROGEN   , ">$OutPath/Intra-Hydrogen_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERHYDROGENSDSD   , ">$OutPath/Inter-Hydrogen_SD_SD_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAHYDROGENSDSD   , ">$OutPath/Intra-Hydrogen_SD_SD_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERHYDROGENSDC   , ">$OutPath/Inter-Hydrogen_SD_C_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAHYDROGENSDC   , ">$OutPath/Intra-Hydrogen_SD_C_$CodePDB.txt") or die "Pb d'ouverture : $!";

              open (INTERHYDROPHOBIC   , ">$OutPath/Inter-Hydrophobic_$CodePDB.txt") or die "Pb d'ouverture : $!";
              open (INTRAHYDROPHOBIC   , ">$OutPath/Intra-Hydrophobic_$CodePDB.txt") or die "Pb d'ouverture : $!";





              #OPEN INPUT CLEANED aligned FILES
              open (FICHIERPDBCLEAN   , "<$PathFileInput") or die "Pb d'ouverture $StructureState/$Subtype/$CodePDB : $!";
              $CountPDB++;   print ".....$flu $StructureState $Subtype $CodePDB processing\n";


              my @LignePDB; my @TdT;

              #Array of array creation
              while ( defined( my $ligne = <FICHIERPDBCLEAN> ) ) {

                @LignePDB = split (" ", $ligne);
                push @TdT, [ @LignePDB ];

              }


                    ## PARALLEL : POTENTIAL INTERACTION COMPUTATIONS ##
                    #Loop : distance calculations --> Put in array all POTENTIAL interactions
                    print "PARALLEL : Search for potential interactions... ";

                    my @PotentialInteracts;
                    my $Parallel_potential = Parallel::ForkManager->new($forks);

                    $Parallel_potential -> run_on_finish ( # called BEFORE the first call to start()
                      sub {
                        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

                          # retrieve data structure from child
                          if (@$data_structure_reference) {  # children are not forced to send anything / Pass Array ref
                            push @PotentialInteracts, $_ foreach @$data_structure_reference;
                          }
                        }
                        );

                        $Parallel_potential->run_on_wait( sub {
                          #print "..."
                        },
                        0.5
                        );


                        my $i; my $j;
                        OUTERLOOP:for ($i=0; $i < $#TdT; $i++)
                        {
                          my @send_to_parent;
                          my $pid = $Parallel_potential->start and next OUTERLOOP;

                            for ($j = $i+1; $j <= $#TdT; $j++) {
                            # 1         	N         	ASP       	A         	11 (CleanNumberAlign)       	-54.024   	-15.147   	-0.774    	HA1       	Proto1    	15 (Number In real pdb file)

                               #DISTANCE CALCULATION IN ANGSTROEM
                               my $Distance = sqrt(($TdT[$i][5] - $TdT[$j][5])**2 + ($TdT[$i][6] - $TdT[$j][6])**2 + ($TdT[$i][7] - $TdT[$j][7])**2);

                                  if ($Distance < 5) {

                                    #Don't consider atomic interaction within a single and same amino acid
                                    unless ( ($TdT[$i][4] == $TdT[$j][4] and $TdT[$i][9] eq $TdT[$j][9]) ) {

                                      my $Potential = "$TdT[$i][1] $TdT[$i][2] $TdT[$i][4] $TdT[$i][8] $TdT[$i][9] old=$TdT[$i][10]   |     $TdT[$j][1] $TdT[$j][2] $TdT[$j][4] $TdT[$j][8] $TdT[$j][9] old=$TdT[$j][10] ->   $Distance";
                                         #  NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

                                      push @send_to_parent, $Potential;
                                    }
                                  }

                                }

                                $Parallel_potential->finish(0, \@send_to_parent);
                              }
                              $Parallel_potential->wait_all_children;

                              print "           --> All potential interactions are retrieved.\n";





                              ################  INTERACTION VALIDATION  ################

                              print "Real interactions screening process...";

                              my %H_AllRealInteractions; #Hash of array : Keys = Interaction Type, Key 2 = Intra or Inter, Value = ref array with all corresponding interactions

                              foreach my $PotentialInteraction (@PotentialInteracts) {

                                #Do we validate the potential interaction ? Explaination more specific distance cutoff + specific atoms in the function
                                my ($validated_interaction, $InteractionType, $IntraOrInter) = calcul_interaction($PotentialInteraction, \%Hash_Interactions_CorrespondingAtoms, \%Hash_Interactions_CorrespondingResidues);

                                if ($validated_interaction) {
                                  push(@{ $H_AllRealInteractions{$InteractionType}{$IntraOrInter} }, $PotentialInteraction);

                                  #ALSO GROUP All Hydrogen and all Repulsive as we have different type of each
                                  if ($InteractionType =~ m/(Hydrogen|Repulsive)/) {
                                    push(@{ $H_AllRealInteractions{$1}{$IntraOrInter} }, $PotentialInteraction);
                                  }
                                }
                              }

                              ################  INTERACTION VALIDATION END ################

                              print "           --> All real interactions have been determined\n";




                              #Print ordered files
                              foreach my $InteractionType (keys %H_AllRealInteractions) {
                                foreach my $intra_inter (keys %{ $H_AllRealInteractions{$InteractionType} } ) {

                                  ## Schwartzian transform to sort the interactions within the hash ##
                                  my @SortedArray = map { $_->[0] }            #retrieve original line
                                  sort {                                     #Sort with Sort::Naturally according to splited line

                                  ncmp($a->[5], $b->[5])   ||
                                  ncmp($a->[4], $b->[4])   ||
                                  $a->[3] <=> $b->[3]    ||    #  Numeric value
                                  ncmp($a->[12], $b->[12]) ||
                                  ncmp($a->[11], $b->[11]) ||
                                  $a->[10] <=> $b->[10]        #  Numeric value

                                }
                                map { [ $_, split " " ] } #Transform each line of array $_->[0] = untransformed data, ->[1] transformed data, -> [1] [2] [3] [4] etc if using split to transform for exemple
                                @{ $H_AllRealInteractions{$InteractionType}{$intra_inter} };
                                #  NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299

                                $H_AllRealInteractions{$InteractionType}{$intra_inter} = \@SortedArray;

                                    foreach my $SortedInteraction ( @SortedArray ) {
                                      #  NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299
                                      my $InfoInteraction = "$intra_inter" . "-" . "$InteractionType";
                                      my $ToPrint = join("\t", $Subtype, $CodePDB, $InfoInteraction, $SortedInteraction);

                                      print_in_files ($ToPrint, $InfoInteraction);
                                    }
                                  }
                                }

                                h_json_and_dumper (\%H_AllRealInteractions, $OutPath, $CodePDB);
                                print "                                                                           have been printed\n";


                                ##### COUNT VALIDATED INTERACTIONS #####
                                my $ref_H_count = count_nb_interactions (\%H_AllRealInteractions);


                                open (my $HandleFileCount, ">$OutPath/CountInteractions_$CodePDB.txt") or die "Pb d'ouverture : $!";
                                open (my $JsonCount,">$OutPath/HashCountInteractions$CodePDB.dat");
                                h_count_json_and_txt (\@AllInteractionTypes, \@Indications, $ref_H_count, $HandleFileCount, $JsonCount);

                                print "                                                                           have been counted\n\n";

                                my $datestring = gmtime();
                                print "__________ $datestring --> $CodePDB DONE / $CountPDB structures DONE __________\n\n\n\n\n";

                }
        }
       }
      }
    }

  print "\n\nProcess DONE for a total of $CountPDB\n\n";












####################################### SUBROUTINES #######################################
####################################### SUBROUTINES #######################################




sub opendir_readdir {
      my $dir = shift;

      opendir(DIR, "$dir") or die "Unable to enter dir $dir:$!\n";
      my @FolderNames = readdir(DIR) or die "Unable to read $dir:$!\n";
      @FolderNames =  grep { ! m/^\./ } @FolderNames;

      return @FolderNames;
                         }





sub h_json_and_dumper {
      my ($H_RealInteractions, $Outpath, $CodePDB) = @_;

      open(HFILE,">$Outpath/HashAllInteractions$CodePDB.dat");
      my $JSON_H = encode_json($H_RealInteractions);   print HFILE $JSON_H;
      close HFILE;

      open(HASHDUMP,">$Outpath/Dumper_HashAllInteractions$CodePDB.txt");
      print HASHDUMP Dumper ($H_RealInteractions);
      close HASHDUMP;
                         }







  sub calcul_interaction {
  my ($Potential_Interaction, $refH_CorrespondingAtoms, $refH_CorrespondingResidues) = @_;

    # input  :  NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299
    # split  :  0   1   2   3     4     5        6      7   8  9   10   11     12    13       14

    ### REFERENCE FOR DISTANCE CUTOFF IN FIRST PAPER ABOUT BIOPHYSIC DIVERGENCE AMONGST HAs



    my @Split = split (" ", $Potential_Interaction);
    foreach (@Split) { $_ =~ s/^\s+|\s+$//g; }

    my $Distancee = $Split[14];
    my $Intra_Inter;

    if ($Split[4] eq $Split[11]) {$Intra_Inter = "Intra";}
      else {$Intra_Inter = "Inter";}


    ###################################### HYDROPHOBIC INTERACTIONS ############################################

    # If atom1 belong to apolar AA's side chain and the same goes for atom2 OR REVERSE ORDER --> Hydrophobic effect
    # Distance should be < 4.2 angstroem


    if ( $Distancee < 4.2 ) {


          if ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomHydrophobic"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduHydrophobic"}} )
          and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomHydrophobic"}} and  any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduHydrophobic"}})) {

     my $Interaction = "Hydrophobic";
     return ($Potential_Interaction, $Interaction, $Intra_Inter);
                      }
    }




    ###################################### HYDROGEN BONDS ############################################

    #Distance hydrogen bonding
    if ($Distancee < 3.5) {


          # We only want sidechain-sidechain or sidechain-backbone hydrogen interactions. (We don't want backbone-backbone)
          unless (( $Split[0] eq "N" and $Split[7] eq "O") or ($Split[0] eq "O" and $Split[7] eq "N")) {

              # We exlude hydrogen interactions present in SALT-BRIDGES
              unless (
                     ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}})
                  and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}}))
                or
                     ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}})
                  and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}}) )
                      ) {

                          # If atom1 is proton donnor and atom2 is proton acceptor OR REVERSE ORDER. (Sidechain-sidechain interaction)
                          if (    ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomDonorProton"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduDonorProton"}})
                               and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomAcceptorProton"}} and  any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduAcceptorProton"}}))
                            or
                                  ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomAcceptorProton"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduAcceptorProton"}})
                               and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomDonorProton"}} and  any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduDonorProton"}}))
                              ) {
                                  my $Interaction = "Hydrogen_SD_SD";
                                  return ($Potential_Interaction, $Interaction, $Intra_Inter);

                                }

                          # (Backbone / Side chain   ou   Side chain / Backbone interactions )
                          # If atom1 is proton donnor (N of backbone) and atom2 is acceptor OR REVERSE ORDER
                          # If atom1 is proton acceptor (O of backbone) and atoms 2 is donor OR REVERSE ORDER
                          elsif (   ($Split[0] eq "N" and  any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomAcceptorProton"}} and  any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduAcceptorProton"}})
                                  or  ($Split[0] eq "O" and  any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomDonorProton"}} and  any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduDonorProton"}})
                                  or  ( any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomAcceptorProton"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduAcceptorProton"}} and $Split[7] eq "N")
                                  or  ( any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomDonorProton"}} and  any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduDonorProton"}} and $Split[7] eq "O")

                                ) {
                                    my $Interaction = "Hydrogen_SD_C";
                                    return ($Potential_Interaction, $Interaction, $Intra_Inter);
                                  }
                                                   }
                                                                                                                                }
                                                }


  ###################################### SALT-BRIDGES ############################################

  #Distance salt-bridges
  if ($Distancee < 4) {

              # If atom1 is positively charged and belong to positively charged AA AND if atom2 is negatively charged and belong to negatively charged AA
              if  ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}})
                and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}}) )
                {
                  my $Interaction = "Ionic";
                  return ($Potential_Interaction, $Interaction, $Intra_Inter);
                }

              # OR REVERSE ORDER
              elsif  ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}})
                      and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}}))
                      {
                        my $Interaction = "Ionic";
                        return ($Potential_Interaction, $Interaction, $Intra_Inter);
                      }
                          }





  ###################################### REPULSIVE INTERACTIONS   ######################################

    #Distance repulsions
    if ($Distancee < 5) {

      # If atom1 is positively charged and belong to positively charged AA AND if atom2 is positively charged and belong to positively charged AA
      if ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}})
          and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicPositive"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicPositive"}}))
          {
            my $Interaction = "RepulsivePositive";
            return ($Potential_Interaction, $Interaction, $Intra_Inter);

          }

       # If atom1 is negatively charged and belong to negatively charged AA AND if atom2 is negatively charged and belong to negatively charged AA
      elsif ((any { $Split[0] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and any { $Split[1] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}})
          and (any { $Split[7] eq $_ } @{$refH_CorrespondingAtoms->{"AtomIonicNegative"}} and any { $Split[8] eq $_ } @{$refH_CorrespondingResidues->{"ResiduIonicNegative"}}))
          {

            my $Interaction = "RepulsiveNegative";
            return ($Potential_Interaction, $Interaction, $Intra_Inter);
                                                                                                                }
          }
  return;
  }







  sub print_in_files {
  my ($ToPrintt, $info) = @_;



   if    ($info eq "Intra-Ionic") {print INTRAIONIC "$ToPrintt\n";}
   elsif ($info eq "Inter-Ionic") {print INTERIONIC "$ToPrintt\n";}

   elsif ($info eq "Intra-Repulsive") {print INTRAREPULSIVE "$ToPrintt\n";}
   elsif ($info eq "Inter-Repulsive") {print INTERREPULSIVE "$ToPrintt\n";}

   elsif ($info eq "Intra-RepulsivePositive") {print INTRAREPULSIVEPOSITIVE "$ToPrintt\n";}
   elsif ($info eq "Inter-RepulsivePositive") {print INTERREPULSIVEPOSITIVE "$ToPrintt\n";}

   elsif ($info eq "Intra-RepulsiveNegative") {print INTRAREPULSIVENEGATIVE "$ToPrintt\n";}
   elsif ($info eq "Inter-RepulsiveNegative") {print INTERREPULSIVENEGATIVE "$ToPrintt\n";}

   elsif ($info eq "Intra-Hydrogen") {print INTRAHYDROGEN "$ToPrintt\n";}
   elsif ($info eq "Inter-Hydrogen") {print INTERHYDROGEN "$ToPrintt\n";}

   elsif ($info eq "Intra-Hydrogen_SD_SD") {print INTRAHYDROGENSDSD "$ToPrintt\n";}
   elsif ($info eq "Inter-Hydrogen_SD_SD") {print INTERHYDROGENSDSD "$ToPrintt\n";}

   elsif ($info eq "Intra-Hydrogen_SD_C") {print INTRAHYDROGENSDC "$ToPrintt\n";}
   elsif ($info eq "Inter-Hydrogen_SD_C") {print INTERHYDROGENSDC "$ToPrintt\n";}

   elsif ($info eq "Intra-Hydrophobic") {print INTRAHYDROPHOBIC "$ToPrintt\n";}
   elsif ($info eq "Inter-Hydrophobic") {print INTERHYDROPHOBIC "$ToPrintt\n";}

   else {print "$info doesn't exist in sub print_in_files --> Need to add it \n";}

   return;
                     }




sub count_nb_interactions {

my $ref_H_AllRealInteractions = shift;
# H = Interaction type -> Intra or Inter -> Array all interacitons
my %H_count;


while ( my ($Type, $h_nest) = each %$ref_H_AllRealInteractions ) {
      foreach my $Intra_Inter (keys %$h_nest) {
            foreach ( @{ $h_nest->{$Intra_Inter} } ) {
                        my $combined_type;
                        my @Split = split (" ", $_);
                     #  NZ LYS 308 HA1 Proto1 old=XX     |     NE2 GLN 62 HA2 Proto1 old=XX  ->   3.31539711045299


                        # EVERY INTERACTION TYPES EXPECT hydrogen and repulsive (cuz we want to count the total hydrogen/repulsive EG Repulsive Positive + Negative)
                        unless ($Type =~ /Hydrogen|Repulsive/) {

                          #Total interaction + Total ITA-p or ITR-p (Intra or Inter)
                          $H_count{"All"}{"Total"}++; $H_count{"All"}{$Intra_Inter}++;
                          #Total interaction type + Total type ITA or OTR
                          $H_count{$Type}{"Total"}++; $H_count{$Type}{$Intra_Inter}++;

                                                                                #All HA1HA1  + All ITA/ITR-p HA1HA1
                              if ($Split[3] eq "HA1" && $Split[10] eq "HA1") { $H_count{"All"}{"HA1-HA1"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA1"}++;
                                                                                #All type HA1HA1 + All type ITA/ITR-p HA1HA1
                                                                               $H_count{$Type}{"HA1-HA1"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                  }

                                                                                #Same for HA1-HA2
                              elsif ( ($Split[3] eq "HA1" && $Split[10] eq "HA2") || ($Split[3] eq "HA2" && $Split[10] eq "HA1") ) {
                                                                                                                                  $H_count{"All"}{"HA1-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                  $H_count{$Type}{"HA1-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                              }

                                                                                #Same for HA2HA2
                              elsif ($Split[3] eq "HA2" && $Split[10] eq "HA2") {    $H_count{"All"}{"HA2-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA2-HA2"}++;
                                                                                     $H_count{$Type}{"HA2-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA2-HA2"}++;
                                                                                  }

                                                              }


                         #SAME PROCESS FOR HYDROGEN (SD-SD - SD-c) / REPULSIVE (Positive - Negative) EXCEPT THAT WE ALSO COUNT THE GLOBAL HYDROGEN and REPULSIVE
                         if ($Type =~ /(Hydrogen_SD_C|Hydrogen_SD_SD|RepulsiveNegative|RepulsivePositive)/) {

                           $combined_type = "Hydrogen" if $Type =~ /Hydrogen_SD_C|Hydrogen_SD_SD/;
                           $combined_type = "Repulsive" if $Type =~ /RepulsiveNegative|RepulsivePositive/;

                           $H_count{"All"}{"Total"}++; $H_count{"All"}{$Intra_Inter}++;
                           $H_count{$Type}{"Total"}++; $H_count{$Type}{$Intra_Inter}++;
                           $H_count{$combined_type}{"Total"}++;    $H_count{$combined_type}{$Intra_Inter}++;

                               if ($Split[3] eq "HA1" && $Split[10] eq "HA1") { $H_count{"All"}{"HA1-HA1"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA1"}++;
                                                                                $H_count{$Type}{"HA1-HA1"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                $H_count{$combined_type}{"HA1-HA1"}++;    $H_count{$combined_type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                   }

                               elsif ( ($Split[3] eq "HA1" && $Split[10] eq "HA2") || ($Split[3] eq "HA2" && $Split[10] eq "HA1") ) {
                                                                                                                                  $H_count{"All"}{"HA1-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                   $H_count{$Type}{"HA1-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                   $H_count{$combined_type}{"HA1-HA2"}++;    $H_count{$combined_type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                               }

                               elsif ($Split[3] eq "HA2" && $Split[10] eq "HA2") {    $H_count{"All"}{"HA2-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA2-HA2"}++;
                                                                                      $H_count{$Type}{"HA2-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA2-HA2"}++;
                                                                                      $H_count{$combined_type}{"HA2-HA2"}++;    $H_count{$combined_type}{$Intra_Inter . "HA2-HA2"}++;
                                                                                   }
                                                                }

                         }
                      }
                    }

  return \%H_count;
                     }





sub h_count_json_and_txt {

my ($ref_ArrayTypes, $ref_ArrayIndications, $ref_HCount, $HandleCountFile, $HandleJson) = @_;



  foreach my $Type (@$ref_ArrayTypes) {
    print $HandleCountFile "\n";

    foreach my $Indication (@$ref_ArrayIndications) {

      my $Info = $Type . "-" . $Indication;

      if ($ref_HCount->{$Type}{$Indication}){
        printf $HandleCountFile ("%-35s", $Info);		# prints 'foo       '
        print $HandleCountFile "\t$ref_HCount->{$Type}{$Indication}\n";
      }

      else {
        printf $HandleCountFile ("%-35s", $Info);		# prints 'foo        '
        print $HandleCountFile "\t0\n";
      }

          }
        }
close $HandleCountFile;


my $JSON_H_count = encode_json($ref_HCount);
print $HandleJson $JSON_H_count;
close $HandleJson;

return;
}
