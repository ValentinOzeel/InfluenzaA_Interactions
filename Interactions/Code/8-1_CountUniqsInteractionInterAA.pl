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




my @AllInteractionTypes = ("All", "Ionic", "Repulsive", "RepulsivePositive", "RepulsiveNegative", "Hydrogen", "Hydrogen_SD_SD", "Hydrogen_SD_C", "Hydrophobic");


my @Indications = ( "Total", "Intra", "Inter", "HA1-HA1", "HA1-HA2", "HA2-HA2",
                      "IntraHA1-HA1", "IntraHA1-HA2", "IntraHA2-HA2",
                      "InterHA1-HA1", "InterHA1-HA2", "InterHA2-HA2" );


my @InteractionTypes = qw(Inter-Ionic             Intra-Ionic
                       Inter-Repulsive            Intra-Repulsive
                       Inter-RepulsivePositive    Intra-RepulsivePositive
                       Inter-RepulsiveNegative    Intra-RepulsiveNegative
                       Inter-Hydrogen             Intra-Hydrogen
                       Inter-Hydrogen_SD_SD       Intra-Hydrogen_SD_SD
                       Inter-Hydrogen_SD_C        Intra-Hydrogen_SD_C
                       Inter-Hydrophobic          Intra-Hydrophobic);




my $CountPDB = 0;



my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");


######################################################## INITIALISATION ############################################################################




my $time_beginning = gmtime();
print "\n\n\n\n__________ $time_beginning : Code start time\n\n\n\n\n\n";





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
              $CountPDB++;
              my $OutPath = "$path_environment/Interactions/Result/Uniqs/$flu/$StructureState/$Subtype/$CodePDB";
              #Create directories for the Output files
              make_path("$OutPath");


                my %H_AllUniqs;

                next if ($CodePDB =~ m/^\./);
                my $PathFileInteractionsPDB = $PathToCodesPDB . "/" . $CodePDB;


              foreach my $interaction_condition (@InteractionTypes) {

                  my $namefile = $interaction_condition . "_" . $CodePDB . ".txt";
                  my $PathINPUT = $PathFileInteractionsPDB . "/" . $namefile;
                  my @split = split("-", $interaction_condition);
                  my $InteractionType = $split[1];
                  my $intra_inter = $split[0];
                  open (INPUT   , "<$PathINPUT") or die "Pb d'ouverture $PathINPUT : $!";


                  my @lines;

                  while ( defined( my $line = <INPUT> ) ) {

                    push(@lines, $line);

                  }


                  my $Uniqs = get_uniq_interaction_interAA(\@lines);


                  $H_AllUniqs{$InteractionType}{$intra_inter} = $Uniqs;


                }



                  h_json_and_dumper (\%H_AllUniqs, $OutPath, $CodePDB);
                  print "                                                                           have been printed\n";


                  ##### COUNT VALIDATED INTERACTIONS #####
                  my $ref_H_count = count_nb_interactions (\%H_AllUniqs);


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




sub get_uniq_interaction_interAA {

     my ($ref_lines) = @_;
     my @lines = @$ref_lines;

     my @uniqs;
    my $count = scalar(@lines);

    while (@lines) {

      my $line = shift(@lines);

       #H5N8	5huf	Intra-Hydrogen_SD_SD	OD1 ASN 23 HA1 Proto1 old=10   |     ND2 ASN 356 HA1 Proto1 old=319 ->   2.86827212795439

       my @split = split(" ", $line);

       my $AA_1 = $split[4]; my $num_1 = $split[5]; my $sub_1 = $split[6]; my $proto_1 = $split[7];
       my $AA_2 = $split[11]; my $num_2 = $split[12]; my $sub_2 = $split[13]; my $proto_2 = $split[14];

       push(@uniqs, $line);

       @lines = grep(!/$AA_1 $num_1 $sub_1 $proto_1.+$AA_2 $num_2 $sub_2 $proto_2/, @lines);
       @lines = grep(!/$AA_2 $num_2 $sub_2 $proto_2.+$AA_1 $num_1 $sub_1 $proto_1/, @lines);

          }


     return \@uniqs;
                        }








sub count_nb_interactions {

my $ref_H_AllUniqs = shift;
# H = Interaction type -> Intra or Inter -> Array all interacitons
my %H_count;



while ( my ($Type, $h_nest) = each %$ref_H_AllUniqs ) {
      foreach my $Intra_Inter (keys %$h_nest) {
            foreach ( @{ $h_nest->{$Intra_Inter} } ) {
                        my $combined_type;
                        my @Split = split (" ", $_);
                     #  H5N1	2fk0	Intra-Ionic	OD1 ASP 193 HA1 Proto1 old=175   |     NZ LYS 256 HA1 Proto1 old=238 ->   3.41614446415839


                        # EVERY INTERACTION TYPES EXPECT hydrogen and repulsive (cuz we want to count the total hydrogen/repulsive EG Repulsive Positive + Negative)
                        unless ($Type =~ /Hydrogen|Repulsive/) {

                          #Total interaction + Total ITA-p or ITR-p (Intra or Inter)
                          $H_count{"All"}{"Total"}++; $H_count{"All"}{$Intra_Inter}++;
                          #Total interaction type + Total type ITA or OTR
                          $H_count{$Type}{"Total"}++; $H_count{$Type}{$Intra_Inter}++;

                                                                                #All HA1HA1  + All ITA/ITR-p HA1HA1
                              if ($Split[6] eq "HA1" && $Split[13] eq "HA1") { $H_count{"All"}{"HA1-HA1"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA1"}++;
                                                                                #All type HA1HA1 + All type ITA/ITR-p HA1HA1
                                                                               $H_count{$Type}{"HA1-HA1"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                  }

                                                                                #Same for HA1-HA2
                              elsif ( ($Split[6] eq "HA1" && $Split[13] eq "HA2") || ($Split[6] eq "HA2" && $Split[13] eq "HA1") ) {
                                                                                                                                  $H_count{"All"}{"HA1-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                  $H_count{$Type}{"HA1-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                              }

                                                                                #Same for HA2HA2
                              elsif ($Split[6] eq "HA2" && $Split[13] eq "HA2") {    $H_count{"All"}{"HA2-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA2-HA2"}++;
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

                               if ($Split[6] eq "HA1" && $Split[13] eq "HA1") { $H_count{"All"}{"HA1-HA1"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA1"}++;
                                                                                $H_count{$Type}{"HA1-HA1"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                $H_count{$combined_type}{"HA1-HA1"}++;    $H_count{$combined_type}{$Intra_Inter . "HA1-HA1"}++;
                                                                                   }

                               elsif ( ($Split[6] eq "HA1" && $Split[13] eq "HA2") || ($Split[6] eq "HA2" && $Split[13] eq "HA1") ) {
                                                                                                                                  $H_count{"All"}{"HA1-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                   $H_count{$Type}{"HA1-HA2"}++; $H_count{$Type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                                   $H_count{$combined_type}{"HA1-HA2"}++;    $H_count{$combined_type}{$Intra_Inter . "HA1-HA2"}++;
                                                                                                                               }

                               elsif ($Split[6] eq "HA2" && $Split[13] eq "HA2") {    $H_count{"All"}{"HA2-HA2"}++; $H_count{"All"}{$Intra_Inter . "HA2-HA2"}++;
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
