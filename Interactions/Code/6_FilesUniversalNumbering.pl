use warnings;
use strict;

use File::Path qw(make_path);
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use JSON;


my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



####### ADD UNIVERSAL NUMERATION (GOT THROUGH 3DCOFFEE ALIGNMENT) TO CLEANED PDB FILES (IMPORTANT FOR CONSERVATION ASSESSMENT BETWEEN DIFFERENT SUBTYPES/CLADES OR EVEN GROUPS)





my %ThreeLetterToOneLetterAACode = (
"ALA" => "A",
"ARG" => "R",
"ASN" => "N",
"ASP" => "D",
"CYS" => "C",
"GLU" => "E",
"GLN" => "Q",
"GLY" => "G",
"HIS" => "H",
"ILE" => "I",
"LEU" => "L",
"LYS" => "K",
"MET" => "M",
"PHE" => "F",
"PRO" => "P",
"SER" => "S",
"THR" => "T",
"TRP" => "W",
"TYR" => "Y",
"VAL" => "V"
);


my $Data_path = "$path_environment/Interactions/Data/";
my $Path_aligment = "Alignments/Alignment_three_proto/";

my $File_handle;

my @Protos = ("Proto1", "Proto2", "Proto3");
my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");





foreach my $flu (@Flu) {

  #OPEN ALIGNMENT FILE AND ACCESS IT
  my $ref_alignment = access_all_aligned_seq($Data_path, $Path_aligment, $flu);
  my $CountPDBdone;

  foreach my $StructureState (@StructureStates) {

   unless ($flu eq "FluB" and $StructureState eq "Intermediates" or $flu eq "FluB" and $StructureState eq "Postfusion"){

    my $path_cleanedsubtypes = "$Data_path/PDB_Clean/Trimeric/$flu/$StructureState";
    my @Subtypes = opendir_readdir ($path_cleanedsubtypes);


      foreach my $Subtype (@Subtypes) {

          next if ($Subtype =~ m/^\./);
          my $path_cleanedfiles = "$path_cleanedsubtypes/$Subtype";
          my @Cleaned_PDB_files = opendir_readdir ($path_cleanedfiles);


          foreach my $Cleaned_PDB_file (@Cleaned_PDB_files) {

            next if ($Cleaned_PDB_file =~ m/^\./);
            next if $Cleaned_PDB_file =~ /^Number/;


            $Cleaned_PDB_file =~ m/PDBClean(.+).txt/;
            my $CodePDB = $1;






            my $Out_printed_file = "$Data_path/PDB_Clean_aligned/Trimeric/$flu/$StructureState/$Subtype";
            make_path($Out_printed_file);

            open (OUTPRINTED, ">", "$Out_printed_file/PDB_Clean_aligned_$CodePDB.txt") or die "Unable to open $Out_printed_file/aligned_$CodePDB.txt :$!\n";
            open (HANDLEPDB, "<", "$path_cleanedfiles/$Cleaned_PDB_file") or die "Unable to enter file $path_cleanedfiles/$Cleaned_PDB_file :$!\n";



            my @PDB_array = ();
            while (<HANDLEPDB>) {
              chomp;
              push @PDB_array, $_;
            }

            my $infos = $flu . " " . $StructureState . " " . $Subtype . " " . $Cleaned_PDB_file . " " . $CodePDB;

            my %H_proto_alignedSeq = ();

                foreach my $Proto (@Protos) {

                  #ACCESS aligned SEQ FOR SPECIFIC PROTOMER
                  my $Seq_aligned_codePDB = get_aligned_sequence ($CodePDB, $Proto, $ref_alignment, $infos);
                  $H_proto_alignedSeq{$Proto} = $Seq_aligned_codePDB;
                }

                 #CREATE HASH WITH KEY : PROTO, Value : Ref array with all corresponding lines
                  my $ref_H_grep_proto = grep_proto (\@PDB_array, \@Protos);

                  #ADD NEW UNIVERSAL NUMERATION
                  my $H_new_numero = alignment_numerotation (\%H_proto_alignedSeq ,\%ThreeLetterToOneLetterAACode, $ref_H_grep_proto, \@Protos, $infos);

                  foreach my $proto (sort(keys %$H_new_numero)) {
                    my $ref_array = $H_new_numero->{$proto};

                      #VERIF IF THE AA IN EACH LINE CORRESPOND TO PREVIOUS PDB CLEAN FILE
                      print_and_verif ($ref_array, $proto, \@PDB_array, $Subtype, $CodePDB);

                      }
                      $CountPDBdone++;


    print "$CodePDB : Done Total : $CountPDBdone Done\n";

    close OUTPRINTED; close HANDLEPDB;




    }
  }
 }
}

  }








##############################################################################################################





    sub opendir_readdir {
          my $dir = shift;

          opendir(DIR, "$dir") or die "Unable to enter dir $dir:$!\n";
          my @FolderNames = readdir(DIR) or die "Unable to read $dir:$!\n";
          @FolderNames =  grep { ! m/^\./ } @FolderNames;
          return @FolderNames;
                             }


  sub access_all_aligned_seq {
         my ($Data_path, $Path_aligment, $flu) = @_;

  local $/ = "\n>";  # read by FASTA record

       #OPEN ALIGNMENT FILE
       open (FILEALIGNMNT, "<", $Data_path.$Path_aligment."Alignment_3proto_Tcoffee_mustang_$flu/final_fasta_mustang_$flu.fa") or die "Unable to enter file $Path_aligment/Alignment_3proto_Tcoffee_mustang_$flu/final_fasta_mustang_$flu.fa :$!\n";
       my $CountPDBdone;

       my @aligned_seq = ();
       while (<FILEALIGNMNT>) {
         chomp;
         push @aligned_seq, $_;
       }
       close FILEALIGNMNT;
       return(\@aligned_seq);


                          }


sub get_aligned_sequence {

  my ($CodePDB, $Proto, $refalignment, $infos) = @_;


  my $aligned;
    foreach (@$refalignment) {

      $aligned = $_ if $_ =~ /$CodePDB$Proto/;  # parse ID as first word in FASTA header
      last if $_ =~ /$CodePDB$Proto/;
    }


      if (defined($aligned)) {
        $aligned =~ s/.{4}Proto\d\n//;  # remove FASTA Alignment ">4k3xProto1\n" header
          $aligned =~ s/\n//g;      # remove endlines
          $aligned =~ s/>//g;      # remove endlines

return $aligned;
}

      if (!defined($aligned)) {
        print "Warning : Can't find $CodePDB$Proto --> $infos $Proto in aligned file\n";}
  }






  sub grep_proto {

  my ($ref_ArrayPDB, $ref_Protos) = @_;

  my %grep_proto;

    foreach my $proto (@$ref_Protos){

          my @grep = grep { $_ =~ /$proto/ } @$ref_ArrayPDB;
          $grep_proto{$proto} = \@grep;
         }

      return (\%grep_proto);
                  }




  sub alignment_numerotation {

  my ($ref_H_proto_alignedSeq, $ref_H_convertAA, $ref_H_grep_proto, $ref_protos, $infos) = @_;


  my %Hash_proto_new_numerotation;

  foreach my $proto (@$ref_protos) {

    my @Array_aligned_pos;

    my $Seq = $ref_H_proto_alignedSeq->{$proto};

    my $ref_arrayPDB_proto = $ref_H_grep_proto->{$proto};


    my $y = 0;


      OUTTER:for (my $i = 0; $i <= length($Seq); $i++ ) {

        #NUMERATION ALIGNMENT
        my $real_pos = $i + 1;

        my $AA_align = substr ($Seq, $i, 1);

        #NEXT IF GAP IN SEQ
        next OUTTER if $AA_align eq "-";


        my $Previous_aa_array; my $Previous_pos_array;

          #WALK TRHOUGH THE LINES OF CLEANED PDB FILE
          INNER:for ($y; $y <= $#$ref_arrayPDB_proto; $y++) {
            my $Line_to_check = $ref_arrayPDB_proto->[$y];
            my @Split = split (" ", $Line_to_check);
            #20        	OE1       	GLN       	A         	6         	350.746   	189.953   	-5.578    	HA1       	Proto1    	6
            #0           1           2           3          4           5              6          7          8           9          10

            foreach (@Split) { $_ =~ s/^\s+|\s+$//g; }
            #REMOVE WHITE SPACE AT THE LEFT AND RIGHT OF A STRING

            my $aa_array = $ref_H_convertAA->{$Split[2]};
            my $pos_array = $Split[4];


                #ENABLE TO BE AWARE OF IF WE CHANGE D  AA AT THIS LOOP ITERATION OR NOT. IF YES WE PASS TO THE NEXT AA OF THE aligned SEQUENCE
                if (!defined($Previous_aa_array) && !defined($Previous_pos_array)) { $Previous_aa_array = $aa_array; $Previous_pos_array = $Split[4]; }
                if ($Previous_aa_array ne $aa_array || $Previous_pos_array != $Split[4]) { $Previous_aa_array = $aa_array; $Previous_pos_array = $Split[4];
                                                                                            last INNER; }


                # WE MODIFY NUMERATION TO THE NEW UNIVERSAL NUMERATION
                if ($AA_align eq $ref_H_convertAA->{$Split[2]}) {


                    my $Brand_new_numerotation = join("\t", $Split[0], $Split[1], $Split[2], $Split[3], $real_pos, $Split[5], $Split[6],
                    $Split[7], $Split[8], $Split[9], $Split[10]);

                    push @Array_aligned_pos, $Brand_new_numerotation;


                }

                else { print "$AA_align $i aligned seq DONT MATCH WITH $aa_array $pos_array : $infos\n"; }

              }


   }

   $Hash_proto_new_numerotation{$proto} = \@Array_aligned_pos;

      }


return \%Hash_proto_new_numerotation;


}





  sub print_and_verif {

  my ($ref_array_to_print, $proto, $ref_arrayPDB, $subtype, $pdbcode) = @_;

  my @Grep_PDB = grep { m/$proto/ } @$ref_arrayPDB;

  my $Count = 0;

    foreach my $line_print (@$ref_array_to_print) {
      my $OKorNot;

          my @Split_to_print = split (" ", $line_print);
          foreach (@Split_to_print) { $_ =~ s/^\s+|\s+$//g; }

          my @Split_PDB = split (" ", $Grep_PDB[$Count]);
          foreach (@Split_to_print) { $_ =~ s/^\s+|\s+$//g; }

          if ($Split_to_print[10] eq $Split_PDB[10] && $Split_to_print[2] eq $Split_PDB[2] && $proto eq $Split_PDB[9]) {

          print OUTPRINTED "$line_print\n";

          }
          else { print "$pdbcode $subtype ISSUE COMPARISON VERIF : $Split_to_print[10] // $Split_PDB[4], $Split_to_print[2] // $Split_PDB[2], $proto // $Split_PDB[9]\n"; }

  $Count++;
  }

  }
