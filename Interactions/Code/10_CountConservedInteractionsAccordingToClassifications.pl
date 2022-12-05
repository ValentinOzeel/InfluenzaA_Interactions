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
my @InteractionTypes = qw(Inter-Ionic             Intra-Ionic
                       Inter-Repulsive            Intra-Repulsive
                       Inter-RepulsivePositive    Intra-RepulsivePositive
                       Inter-RepulsiveNegative    Intra-RepulsiveNegative
                       Inter-Hydrogen             Intra-Hydrogen
                       Inter-Hydrogen_SD_SD       Intra-Hydrogen_SD_SD
                       Inter-Hydrogen_SD_C        Intra-Hydrogen_SD_C
                       Inter-Hydrophobic          Intra-Hydrophobic);



my @AllInteractionTypes = ("All", "Ionic", "Repulsive", "RepulsivePositive", "RepulsiveNegative", "Hydrogen", "Hydrogen_SD_SD", "Hydrogen_SD_C", "Hydrophobic");


my @Indications = ( "Total", "Intra", "Inter", "HA1-HA1", "HA1-HA2", "HA2-HA2",
                       "IntraHA1-HA1", "IntraHA1-HA2", "IntraHA2-HA2",
                       "InterHA1-HA1", "InterHA1-HA2", "InterHA2-HA2" );


##### COUNT CONSERVED INTERACTIONS #####
########## %H_count_conservedhighest HASH IS THE ONLY ONE TO CONSIDER (>= 90% CONSERVATION).
########## %H_count_conservedall hash take into account lower conservation rates.


my %H_count_conservedall;
my %H_count_conservedhighest;
my %H_variability_possibilities;

foreach my $flu (@Flu) {
print $flu, "\n";
 foreach my $state (@StructureStates) {
  unless ($flu eq "FluB" and $state eq "Intermediates" or $flu eq "FluB" and $state eq "Postfusion"){
  print $state, "\n";
   foreach my $type (@InteractionTypes) {

   my $In_P = "$path_environment/Interactions/Result/Conservation/$flu/$state/$type/JSON/Conservation_values_$type.dat";

   if (-e $In_P) {

   my $H_conserv = get_json_data_hardsaved ($In_P);

   my $ref_H_conserv = dclone($H_conserv);



        my @interactions_keys = (keys %{ $ref_H_conserv->{$flu}{$state}{$type} } );


########## CALCUL CONSERVED INTERACTION NB COMBINATIONS

        my %H_typeconserv;

        foreach (@interactions_keys) {
        #WE GROUP DIFFERENT CONSERVED INTERACTIONS WHICH ARE RELATED (SAME POSITIONS BUT MUTATION PRESERVING TYPE) = ALL POSSIBLE AA COMBINATIONS PRESERVING INTERACTION TYPE
        #PREVENT RUNNING MULTIPLE TIME INTERACTION BETWEEN TWO IDENTICAL POSITIONS (but perhaps different amino acid)

          $H_typeconserv{$_} = [];

          #Inter-Ionic ARG 340 HA1 ASP 454 HA2
          my @Split = split (" ", $_);
          foreach (@Split) { $_ =~ s/^ +| +$//g; }

          my $regex1 = "$Split[2] HA.+ $Split[5] HA";
          my $regex2 = "$Split[5] HA.+ $Split[2] HA";

          my $index = 0;
          $index++ until $interactions_keys[$index] eq $_;
          splice(@interactions_keys, $index, 1);

          my @grep = grep { /$regex1|$regex2/ } @interactions_keys;

                foreach my $Grep (@grep) {
                    my @Splt = split (" ", $Grep);

                    if (($Split[2] == $Splt[5] && $Split[5] == $Splt[2])|| ($Split[2] == $Splt[2] && $Split[5] == $Splt[5])) {

                                          delete( $ref_H_conserv->{$flu}{$state}{$type}{$Grep} );

                                          push @{ $H_typeconserv{$_} }, $Grep;
                                        }
                                  }
        }


        #COUNT NUMBER OF CONSERVED INTERACTION ACCROSS CLASSIFICATIONS + TOTAL NB OF AA COMBINATIONS RETRIEVED (Enabling to preserve interaction type)
        while (my ($interaction, $correspondings_interacts) = each (%H_typeconserv) ) {

          my $number_possibilities = 1 + scalar( @{ $correspondings_interacts } );

          #COUNT DETAILED NUMBER OF COMBINATIONS (POSSIBILITIES) ACCORDING TO DIFFERENT CONDITIONS (HA1HA1, HA1HA2 , ITR-p ...) and add results in H_variability_possibilities
          count_nb_interactions_conserv ($flu, $state, undef, $type, $interaction, undef, undef, \%H_variability_possibilities, undef, $number_possibilities);
              }

########## END CALCUL CONSERVED INTERACTION NB COMBINATIONS





        count_conserved ($ref_H_conserv, $flu, $state, $type, \%H_count_conservedall, \%H_count_conservedhighest);




          }

    }



    foreach my $type (@AllInteractionTypes) {
        foreach my $indication (@Indications) {



          my $State_data_array = "Prefusion";    $State_data_array = "Postfusion_Intermediary" if $state =~ m/Intermediates|Postfusion/;
          my $Path_class = "$path_environment/Interactions/Data/ArrayTableauPDB/$flu/$State_data_array/HashAllConditionsGrep_$State_data_array.dat";
          my $ref_H_category = get_json_data_hardsaved ($Path_class);
          delete $ref_H_category->{"Subtypes_Hosts"} if exists $ref_H_category->{"Subtypes_Hosts"};
          @{ $ref_H_category->{"Globally"}{"Globally"} } = [];

          my $copy_H_category = dclone($ref_H_category);

          while (keys(%$copy_H_category)) {

            my ($last_key, $ref_array_v, $ref_key_list) = hash_walk_through ($copy_H_category, [], \&deletekey);
            my $first_key = $ref_key_list->[0]; my $second_key = $ref_key_list->[1];



            if (defined($first_key) && defined($second_key)) {


                    if (exists $H_count_conservedall{$flu}{$state}{$type}{$indication}{$first_key}{$second_key}) { }
                    else { $H_count_conservedall{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} = 0; }


                    if (exists $H_count_conservedhighest{$flu}{$state}{$type}{$indication}{$first_key}{$second_key}) { }
                    else { $H_count_conservedhighest{$flu}{$state}{$type}{$indication}{$first_key}{$second_key} = 0; }

                  }
              }

          }
        }


    }
  }
}





#print Dumper \%H_variability_possibilities;
#print Dumper \%H_count_conservedhighest;











##### CREATE FILES (PRISM FORMAT) TO PERFORM FURTHER ANALYSES #####
##### NB CONSERVED INTERACTIONS ACCROSS CLASSIFICATIONS
print "\n\n    Creating files ...\n\n";
my $index = 0;

my @Hashes_ref = (\%H_count_conservedall, \%H_count_conservedhighest);
my @infos_folder = ("Count_conserved", "Count_highest_conserved");


for (@Hashes_ref) {


my $info_folder = $infos_folder[$index];
print $info_folder, "\n";
$index++;

foreach my $flu (@Flu) {
  print " $flu\n";
  foreach my $state (@StructureStates) {
    unless ($flu eq "FluB" and $state eq "Intermediates" or $flu eq "FluB" and $state eq "Postfusion"){
    print "  $state\n";

    foreach my $type (nsort keys %{ $_->{$flu}{$state} }) {
    print "    $type\n";
      foreach my $indication (nsort keys %{ $_->{$flu}{$state}{$type} }) {
      my $Type_indic = $type . " " . $indication;
        print "      $indication\n";

          foreach my $first_key (nsort keys %{ $_->{$flu}{$state}{$type}{$indication} }) {
            my %second_keys_values = ();

                    foreach my $second_key (nsort keys %{ $_->{$flu}{$state}{$type}{$indication}{$first_key} }) {

                      my $value_conserv = $_->{$flu}{$state}{$type}{$indication}{$first_key}{$second_key};

                      $second_keys_values{$second_key} = $value_conserv;
                      }

                          my $namefile = $Type_indic . " " . $first_key . ".txt";
                          my $path_out = "$path_environment/Interactions/Result/$info_folder/$flu/$state/$first_key/$Type_indic";
                          make_path($path_out);
                          open (FILETOPRINT   , ">$path_out/$namefile") or die "Pb d'ouverture $path_out/$namefile : $!";


#                            PRINT HEADER`
                          printf FILETOPRINT ("%-15s", $_) foreach (nsort keys %second_keys_values);
                          print FILETOPRINT "\n";

                          print "        PROCESSING : ", $first_key, "\n";





                          my $number_second_key = scalar (keys %second_keys_values);
                          my $nb_keydone = 0;


                          WHILE:while ($nb_keydone != $number_second_key) {


                            foreach my $key (nsort keys %second_keys_values) {


                              my $to_print = $second_keys_values{$key};
                              printf FILETOPRINT ("%-15s", $to_print);

                              $nb_keydone++;


                          }
                            printf FILETOPRINT "\n";
                        }

      }
    print "\n\n";
     }
    }
    }
   }
  }
}






##### CREATE FILES (PRISM FORMAT) TO PERFORM FURTHER ANALYSES #####
##### NB COMBINATIONS POSSIBILITIES
foreach my $flu (@Flu) {
  foreach my $state (@StructureStates) {
    unless ($flu eq "FluB" and $state eq "Intermediates" or $flu eq "FluB" and $state eq "Postfusion"){
    foreach my $type (nsort keys %{ $H_variability_possibilities{$flu}{$state} }) {
      foreach my $indication (nsort keys %{ $H_variability_possibilities{$flu}{$state}{$type} }) {
      my $Type_indic = $type . " " . $indication;
      my %H_for_file;

              foreach my $key_nb (nsort keys %{ $H_variability_possibilities{$flu}{$state}{$type}{$indication} }) {
                  $H_for_file{$key_nb} = $H_variability_possibilities{$flu}{$state}{$type}{$indication}{$key_nb};
              }

        my $namefile = $Type_indic . " " . "Possibilities" . ".txt";
        my $path_out = "$path_environment/Interactions/Result/Variability_Possibilities/$flu/$state/$Type_indic";
        make_path($path_out);
        open (FILETOPRINT   , ">$path_out/$namefile") or die "Pb d'ouverture $path_out/$namefile : $!";

  #      PRINT HEADER`
        printf FILETOPRINT ("%-15s", $_) foreach (nsort keys %H_for_file);
        print FILETOPRINT "\n";

        my $number_key = scalar (keys %H_for_file);
        my $nb_keydone = 0;

        WHILE:while ($nb_keydone != $number_key) {
          foreach my $key (nsort keys %H_for_file) {

            my $to_print = $H_for_file{$key};
            printf FILETOPRINT ("%-15s", $to_print);
            $nb_keydone++;
        }
          printf FILETOPRINT "\n";
      }

        }
      }
     }
    }
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
        # Keep track of the hierarchy of keys, in case
        # our callback needs it.
        push @$key_list, $key;
        if ($key_list->[0] ne $key) {
        if ( (any {$key_list->[0] eq $_} (@keys)) && (any {$key eq $_} (@keys)) ) { shift @$key_list;  }
            }

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






sub count_conserved {

my ($ref_H_conserv, $flu, $state, $type, $ref_H_countconserv, $ref_H_countconservhighest) = @_;


   foreach my $interaction (keys %{ $ref_H_conserv->{$flu}{$state}{$type} } ) {
        foreach my $conserv_info (keys %{ $ref_H_conserv->{$flu}{$state}{$type}{$interaction} } ) {
            foreach my $key1 (keys %{ $ref_H_conserv->{$flu}{$state}{$type}{$interaction}{$conserv_info} } ) {
                foreach my $key2 (keys %{ $ref_H_conserv->{$flu}{$state}{$type}{$interaction}{$conserv_info}{$key1} } ) {


               count_nb_interactions_conserv ($flu, $state, $conserv_info, $type, $interaction, $key1, $key2, $ref_H_countconserv, $ref_H_countconservhighest, );

             }
           }
         }
       }



return;

}




sub count_nb_interactions_conserv {

my ($flu, $state, $conserv_info, $interact_type, $uniq_interact, $first_key, $second_key, $ref_H_count_conserved, $ref_H_count_conservedhighest, $number_possibilities) = @_;
# H = Interaction type -> Intra or Inter -> Array all interacitons

my @Split   = split (" ", $uniq_interact);
#  Intra-Ionic LYS 308 HA1 GLU 316 HA1

my $Split_3 = $Split[3];
my $Split_6 = $Split[6];

my @Split_type = split ("-", $interact_type);
my $Intra_Inter = $Split_type[0];
my $Type = $Split_type[1];


   #USE THIS SUB TO COUNT POSSIBILITIES HERE $ref_H_count_conserved IS IN FACT \%H_variability_possibilities
   #OR 2ND USE CASE IS COUNT CONSERVED INTERACTIONS Conserved + Highly conserved
   detail_count_interaction ($ref_H_count_conserved, $flu, $state, $first_key, $second_key, $Split_3, $Split_6, $Intra_Inter, $Type, $number_possibilities);




   return if !defined($conserv_info);

   my @highest_conserv = qw(mostly_conserved fully_conserved);

   if ( grep { m/$conserv_info/ } @highest_conserv ) {   detail_count_interaction ($ref_H_count_conservedhighest, $flu, $state, $first_key, $second_key, $Split_3, $Split_6, $Intra_Inter, $Type, );   }



  return;
                         }









sub detail_count_interaction {

my ($ref_H, $flu, $state, $first_key, $second_key, $Split_3, $Split_6, $Intra_Inter, $Type, $number_possibilities) = @_;



if (!defined($first_key) && !defined($second_key)) {


  $ref_H->{$flu}{$state}{"All"}{"Total"}{"nb:".$number_possibilities}++;
  $ref_H->{$flu}{$state}{"All"}{$Intra_Inter}{"nb:".$number_possibilities}++;

  $ref_H->{$flu}{$state}{$Type}{"Total"}{"nb:".$number_possibilities}++;
  $ref_H->{$flu}{$state}{$Type}{$Intra_Inter}{"nb:".$number_possibilities}++;


      if ($Split_3 eq "HA1" && $Split_6 eq "HA1") { $ref_H->{$flu}{$state}{"All"}{"HA1-HA1"}{"nb:".$number_possibilities}++;
                                                      $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA1-HA1"}{"nb:".$number_possibilities}++;

                                                       $ref_H->{$flu}{$state}{$Type}{"HA1-HA1"}{"nb:".$number_possibilities}++;
                                                       $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA1-HA1"}{"nb:".$number_possibilities}++;
                                                          }

      elsif ( ($Split_3 eq "HA1" && $Split_6 eq "HA2") || ($Split_3 eq "HA2" && $Split_6 eq "HA1") ) {
                                                          $ref_H->{$flu}{$state}{"All"}{"HA1-HA2"}{"nb:".$number_possibilities}++;
                                                          $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA1-HA2"}{"nb:".$number_possibilities}++;

                                                          $ref_H->{$flu}{$state}{$Type}{"HA1-HA2"}{"nb:".$number_possibilities}++;
                                                          $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA1-HA2"}{"nb:".$number_possibilities}++;
                                                                                                      }

      elsif ($Split_3 eq "HA2" && $Split_6 eq "HA2") {  $ref_H->{$flu}{$state}{"All"}{"HA2-HA2"}{"nb:".$number_possibilities}++;
                                                            $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA2-HA2"}{"nb:".$number_possibilities}++;

                                                             $ref_H->{$flu}{$state}{$Type}{"HA2-HA2"}{"nb:".$number_possibilities}++;
                                                             $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA2-HA2"}{"nb:".$number_possibilities}++;
                                                          }




}






else {


$ref_H->{$flu}{$state}{"All"}{"Total"}{$first_key}{$second_key}++;
$ref_H->{$flu}{$state}{"All"}{$Intra_Inter}{$first_key}{$second_key}++;

$ref_H->{$flu}{$state}{$Type}{"Total"}{$first_key}{$second_key}++;
$ref_H->{$flu}{$state}{$Type}{$Intra_Inter}{$first_key}{$second_key}++;


    if ($Split_3 eq "HA1" && $Split_6 eq "HA1") { $ref_H->{$flu}{$state}{"All"}{"HA1-HA1"}{$first_key}{$second_key}++;
                                                    $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA1-HA1"}{$first_key}{$second_key}++;

                                                     $ref_H->{$flu}{$state}{$Type}{"HA1-HA1"}{$first_key}{$second_key}++;
                                                     $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA1-HA1"}{$first_key}{$second_key}++;
                                                        }

    elsif ( ($Split_3 eq "HA1" && $Split_6 eq "HA2") || ($Split_3 eq "HA2" && $Split_6 eq "HA1") ) {
                                                        $ref_H->{$flu}{$state}{"All"}{"HA1-HA2"}{$first_key}{$second_key}++;
                                                        $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA1-HA2"}{$first_key}{$second_key}++;

                                                        $ref_H->{$flu}{$state}{$Type}{"HA1-HA2"}{$first_key}{$second_key}++;
                                                        $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA1-HA2"}{$first_key}{$second_key}++;
                                                                                                    }

    elsif ($Split_3 eq "HA2" && $Split_6 eq "HA2") {  $ref_H->{$flu}{$state}{"All"}{"HA2-HA2"}{$first_key}{$second_key}++;
                                                          $ref_H->{$flu}{$state}{"All"}{$Intra_Inter . "HA2-HA2"}{$first_key}{$second_key}++;

                                                           $ref_H->{$flu}{$state}{$Type}{"HA2-HA2"}{$first_key}{$second_key}++;
                                                           $ref_H->{$flu}{$state}{$Type}{$Intra_Inter . "HA2-HA2"}{$first_key}{$second_key}++;
                                                        }

      }




return;

      }






      sub h_json_and_dumper {
            my ($H_tosave, $Outpath, $Info) = @_;

            my $path_json = "$Outpath/JSON"; my $path_dump = "$Outpath/DUMPER";
            make_path($path_json);           make_path($path_dump);

            open(HFILE,">$path_json/$Info.dat") or die "Impossible to open : $path_json/$Info.dat : $!";;
            my $path_json_file = "$path_json/$Info.dat";
            my $JSON_H = encode_json($H_tosave);
            print HFILE $JSON_H;
            close HFILE;

            print "\n------> Hash $Info : Saved\n";

            open(HASHDUMP,">$path_dump/$Info.txt");
            print HASHDUMP Dumper ($H_tosave);
            close HASHDUMP;

        return $path_json_file;
                               }
