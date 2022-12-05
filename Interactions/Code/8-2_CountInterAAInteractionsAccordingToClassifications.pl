use warnings;
use strict;

use List::Util qw(sum uniq any max reduce);

use JSON;

use Sort::Naturally;

use Storable qw(dclone);

use File::Path qw(make_path);

use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Parallel::ForkManager;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;




my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment




my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");


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




my @list_subtypes = ();
my @List_all_codePDB = ();

my @conservation_info = qw(fairly_conserved decently_conserved highly_conserved mostly_conserved fully_conserved);

my $info = Sys::Info->new;
my $cpu  = $info->device( 'CPU' );
my $NumberCPU = $cpu->count;


############### GET LIST ALL SUBTYPE AND CODE PDB ###############

            foreach my $flu (@Flu) {
              foreach my $StructureState (@StructureStates){
                unless ($flu eq "FluB" and $StructureState eq "Intermediates" or $flu eq "FluB" and $StructureState eq "Postfusion"){
                my $dirInputSubtype = "$path_environment/Interactions/Result/Uniqs/$flu/$StructureState";
                my @SubtypesInputToProcess = opendir_readdir ($dirInputSubtype);

                    foreach my $Subtype (@SubtypesInputToProcess) {

                    next if ($Subtype =~ m/^\./);
                    push @list_subtypes, $Subtype;
                    my $PathToCodesPDB = $dirInputSubtype . "/" . $Subtype;
                    my @CodesPDBInputToProcess = opendir_readdir ($PathToCodesPDB);

                        foreach my $CodePDB (@CodesPDBInputToProcess) {
                        next if ($CodePDB =~ m/^\./);
                        push @List_all_codePDB, $CodePDB;

                    }
                  }
                 }
                }
              }
      @list_subtypes = uniq @list_subtypes;
      @List_all_codePDB = uniq @List_all_codePDB;

############### GET LIST ALL SUBTYPE AND CODE PDB ###############







##### CREATE HASH COUNT IN FUNCTION OF THE DIFFERENT CONDITONS : GROUP, CLADE, SUBTYPE....  #####
print "\n\nCalculations : hash count ...\n\n";
my %H_count_normal;
my %H_proportionvstotal_normal;
my %H_proportionvstotaltype_normal;
foreach my $flu (@Flu) {
  foreach my $state (@StructureStates) {
  unless ($flu eq "FluB" and $state eq "Intermediates" or $flu eq "FluB" and $state eq "Postfusion"){
  my $State_data_array;

    foreach my $type (@AllInteractionTypes) {
      foreach my $indication (@Indications) {


    my $State_data_array = "Prefusion";    $State_data_array = "Postfusion_Intermediary" if $state =~ m/Intermediates|Postfusion/;
    my $Path_class = "$path_environment/Interactions/Data/ArrayTableauPDB/$flu/$State_data_array/HashAllConditionsGrep_$State_data_array.dat";
    my $ref_H_category = get_json_data_hardsaved ($Path_class);
    delete $ref_H_category->{"Subtypes_Hosts"} if exists $ref_H_category->{"Subtypes_Hosts"};

    my $copy_H_category = dclone($ref_H_category);


    while (keys(%$copy_H_category)) {
      my ($last_key, $ref_array_v, $ref_key_list) = hash_walk_through ($copy_H_category, [], \&deletekey);

      my $first_key = $ref_key_list->[0]; my $second_key = $ref_key_list->[1];


      LOOP:foreach $_ (@$ref_array_v) {

          my @split = split (/\|/, $_);

          my $PDB_code = lc $split[0];
          my $subtype = $split[1];
          $subtype =~ s/B_//g if $flu eq "FluB";

          my $In_path = "$path_environment/Interactions/Result/Uniqs/$flu/$state/$subtype/$PDB_code/HashCountInteractions$PDB_code.dat";



            if (-e $In_path) {

                    my $ref_H_count = get_json_data_hardsaved ($In_path);

                    my $value = $ref_H_count->{$type}{$indication};
                    my $value_all_interactions = $ref_H_count->{"All"}{"Total"};
                    my $value_all_type = $ref_H_count->{$type}{"Total"};





                    if (!defined $value) {
                                            push @{ $H_count_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, 0;
                                            push @{ $H_proportionvstotal_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, 0;
                                            push @{ $H_proportionvstotaltype_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, 0;
                                            }

                    else {
                          my $value_proportion = $value/$value_all_interactions;
                          my $value_proportiontype = $value/$value_all_type;



                          push @{ $H_count_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, $value;


                          my $proportion_vstotal = sprintf("%.5f", $value_proportion);
                          push @{ $H_proportionvstotal_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, $proportion_vstotal;


                          my $proportion_vstypetotal = sprintf("%.5f", $value_proportiontype);
                          push @{ $H_proportionvstotaltype_normal{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} }, $proportion_vstypetotal;

                            }

                    }
      }
     }

        }
      }
     }
    }
  }


  my $HandleJson_normal = "$path_environment/Interactions/Result/Uniqs/Count_normal";
  my $HandleJson_proportionvstotal = "$path_environment/Interactions/Result/Uniqs/Count_proportionvstotal";
  my $HandleJson_proportionvstotaltype = "$path_environment/Interactions/Result/Uniqs/Count_proportionvstotaltype";

  save_json(\%H_count_normal, $HandleJson_normal, "Count_H_normal");
  save_json(\%H_proportionvstotal_normal, $HandleJson_proportionvstotal, "Count_H_proportionvstotal_normal");
  save_json(\%H_proportionvstotaltype_normal, $HandleJson_proportionvstotaltype, "Count_H_proportionvstotaltype_normal");


  print "\n\n_________________________________________\n\nHASH COUNT NORMAL / PROPORTION AND PROPORTIONTYPE HAVE BEEN SUCCESSFULLY SAVED\n\n_________________________________________\n\n";

##### END CREATE HASH COUNT IN FUNCTION OF THE DIFFERENT CONDITONS : GROUP, CLADE, SUBTYPE.... #####





my @Hashes_ref = (\%H_count_normal, \%H_proportionvstotal_normal, \%H_proportionvstotaltype_normal);
my @infos_folder = ("Count_normal", "Count_proportionvstotal", "Count_proportionvstotaltype");

##### CREATE COUNT FILES (PRISM FORMAT) TO PERFORM STATISTICS #####
print "  Creating files ...\n\n";
my $index = 0;

for (@Hashes_ref) {

my $info_folder = $infos_folder[$index];
print $info_folder, "\n";
$index++;

foreach my $flu (@Flu) {
  print " $flu\n";
  foreach my $state (@StructureStates) {
    unless ($flu eq "FluB" and $state eq "Intermediates" or $flu eq "FluB" and $state eq "Postfusion"){  
    print "  $state\n";

    foreach my $type (keys %{ $_->{$flu}{$state} }) {
    print "    $type\n";
      foreach my $indication (keys %{ $_->{$flu}{$state}{$type} }) {
      my $Type_indic = $type . " " . $indication;
        print "      $indication\n";

          foreach my $first_key (keys %{ $_->{$flu}{$state}{$type}{$indication} }) {
            my %second_keys_values = ();

                    foreach my $second_key (keys %{ $_->{$flu}{$state}{$type}{$indication}{$first_key} }) {
                      my @array = ();
                      push @array, $_ foreach @{ $_->{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} };

                      $second_keys_values{$second_key} = \@array;
                      }

                          my $namefile = $Type_indic . " " . $first_key . ".txt";
                          my $path_out = "$path_environment/Interactions/Result/Uniqs/$info_folder/$flu/$state/$first_key/$Type_indic";
                          make_path($path_out);
                          open (FILETOPRINT   , ">$path_out/$namefile") or die "Pb d'ouverture $path_out/$namefile : $!";



                          my $number_second_key = scalar (keys %second_keys_values);
                          my @key_done = ();
                          my $nb_keydone = scalar(@key_done);


#                            PRINT HEADER
                          printf FILETOPRINT ("%-15s", $_) foreach (nsort keys %second_keys_values);
                          print FILETOPRINT "\n";

                          print "        PROCESSING : ", $first_key, "\n";


                          WHILE:while ($nb_keydone != $number_second_key) {
                            my $counter = 0;

                            foreach my $key (nsort keys %second_keys_values) {
                              my $nb_values = scalar(@{$second_keys_values{$key}});
                              $counter++;

                              if ( $nb_values > 0 ) { my $to_print = shift @{$second_keys_values{$key}};
                                                                   printf FILETOPRINT ("%-15s", $to_print);
                                                              }
                                                elsif ($nb_values == 0 ) {
                                                                        push @key_done, $key unless any{$key eq $_} @key_done;
                                                                        my $nb = scalar(@key_done);
                                                                        printf FILETOPRINT ("%-15s", "-");
                                                                        last WHILE if $nb == $number_second_key;
                                                                      }

                              if ($counter == $number_second_key) {
                                          print FILETOPRINT "\n";
                                          next WHILE;
                                  }

                          }
                        }

      }
    print "\n\n";
     }
    }
    }
   }
  }


}

  print "\n\n_________________________________________\n\nFILES FOR STATISTICAL ANALYSIS HAVE BEEN SUCCESSFULLY CREATED\n\n_________________________________________\n\n";
##### END CREATE FILES TO PERFORM STATISTICS #####
















###################################### SUBS ################################################















       sub save_json {
             my ($H_ref, $Outpath, $Info) = @_;

             open(HFILE,">$Outpath/$Info.dat");
             my $JSON_H = encode_json($H_ref);
             print HFILE $JSON_H;
             close HFILE;



         return;
                                }



       sub opendir_readdir {
             my $dir = shift;

             opendir(DIR, "$dir") or die "Unable to enter dir $dir:$!\n";
             my @FolderNames = readdir(DIR) or die "Unable to read $dir:$!\n";
             @FolderNames =  grep { ! m/^\./ } @FolderNames;
             return @FolderNames;
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

    my @keys = qw(SubtypesSameHA Subtypes Hosts HA_Clades HA_Groups);
    my $KEY; my $VALUE; my $ref_KEY_LIST;

    while (my ($key, $value) = each %$hash) {
        # Keep track of the hierarchy of keys, in case
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




sub count_nb_interactions {

my ($flu, $state, $conserv_info, $interact_type, $uniq_interact, $first_key, $second_key, $ref_H_count_conserved) = @_;
# H = Interaction type -> Intra or Inter -> Array all interacitons

my @Split   = split (" ", $uniq_interact);
#  Intra-Ionic LYS 308 HA1 GLU 316 HA1
my @Split_type = split ("-", $interact_type);
my $Intra_Inter = $Split_type[0];
my $Type        = $Split_type[1];
my $conserv;

if    ($conserv_info eq "fairly_conserved" || $conserv_info eq "decently_conserved") {$conserv = "f_d_conserved";}
elsif ($conserv_info eq "highly_conserved" || $conserv_info eq "mostly_conserved") {$conserv = "h_m_conserved";}
elsif ($conserv_info eq "fully_conserved") {$conserv = "fully_conserved";}



  $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{"Total"}{$first_key}{$second_key}++;
  $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{$Intra_Inter}{$first_key}{$second_key}++;

  $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{"Total"}{$first_key}{$second_key}++;
  $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{$Intra_Inter}{$first_key}{$second_key}++;


      if ($Split[3] eq "HA1" && $Split[6] eq "HA1") { $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{"HA1-HA1"}{$first_key}{$second_key}++;
                                                      $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{$Intra_Inter . "HA1-HA1"}{$first_key}{$second_key}++;

                                                       $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{"HA1-HA1"}{$first_key}{$second_key}++;
                                                       $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{$Intra_Inter . "HA1-HA1"}{$first_key}{$second_key}++;
                                                          }

      elsif ( ($Split[3] eq "HA1" && $Split[6] eq "HA2") || ($Split[3] eq "HA2" && $Split[6] eq "HA1") ) {
                                                          $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{"HA1-HA2"}{$first_key}{$second_key}++;
                                                          $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{$Intra_Inter . "HA1-HA2"}{$first_key}{$second_key}++;

                                                          $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{"HA1-HA2"}{$first_key}{$second_key}++;
                                                          $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{$Intra_Inter . "HA1-HA2"}{$first_key}{$second_key}++;
                                                                                                      }

      elsif ($Split[3] eq "HA2" && $Split[6] eq "HA2") {  $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{"HA2-HA2"}{$first_key}{$second_key}++;
                                                            $ref_H_count_conserved->{$flu}{$state}{$conserv}{"All"}{$Intra_Inter . "HA2-HA2"}{$first_key}{$second_key}++;

                                                             $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{"HA2-HA2"}{$first_key}{$second_key}++;
                                                             $ref_H_count_conserved->{$flu}{$state}{$conserv}{$Type}{$Intra_Inter . "HA2-HA2"}{$first_key}{$second_key}++;
                                                          }
  return;
                         }
