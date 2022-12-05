use warnings;
use strict;

use File::Path qw( make_path );


############# CODE PDB CLEAN for structures which DON'T delimit HA1 (Residues 1 to ~325) and HA2 (Residues 1 to ~150).
############# THEY DO ONLY PRESENT HA0 NUMEROTATION so NEED TO IDENTIFY HA1 AND HA2 RESIDUES

my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment




##########################


# Hash of : identified chain names combinations (KEYS) within raw PDB files => corresponding PDB files (Value = ref list)
## NUMERO OF FIRST HA2 NUMBER IS INDICATED AFTER PDB CODE by "_XXX" identified according to (next line)
## CLEAVAGE HA0 :  HA1...R/GLFGAIA....HA2   The first amino acid after "/" (representing cleavage site) is the first aa of HA2

my %HashPDBinfo = (
"_A_ H1N1 A/D/F" => [ "6d8w_330" ],

"_A_ H1N1 A/B/C" => [ "6n41_330", "6ona_330", "6osr_330" ],

"_A_ H2N1 A/B/C" => [ "2wr5_330", "2wr3_330", "2wr4_330" ],

"_A_ H2N2 A/B/C" => [ "2wr7_330", "2wrd_330", "2wre_330", "2wrf_330", "2wr0_330", "2wr2_330", "2wr1_330", "2wrc_330", "2wrb_330" ],

"_A_ H3N2 A/B/C" => [ "3whe_330", "5kaq_330", "6mxu_330", "6mym_330", "6mzk_330", "6n08_330", "6pdx_330", "6wxb_330" ],

"_A_ H3N2 A/C/D" => [ "6p6p_330" ],

"_A_ H3N2 A/B C/D E/F" => [ "5kuy_330" ],

"_A_ H3N8 B/C/D" => [ "4wa1_330", "4wa2_330" ],

"_A_ H4N6 A/B/C" => [ "5xl2_328", "5xl8_328", "5xl9_328" ],

"_A_ H4N6 A/D B/E C/F" => [ "5y2m_328" ],

"_A_ H5N1 A/B C/D E/F" => [ "4kdm_335", "4kdn_335", "4kdo_335" ],

"_A_ H5N2 A/B/C" => [ "5ykc_330" ],

"_A_ H6N1 A/B C/D E/F" => [ "4yy0_330" ],

"_A_ H6N2 A/B/C" => [ "4wsr_329", "4wss_329" ],

"_A_ H7N9 A/B C/D E/F" => [ "4lcx_322", "6ii8_322", "6id2_322" ],

"_A_ H7N9 A/E/G" => [ "6idd_322" ],

"_A_ H10N8 A/B C/D E/F" => [ "4qy0_324", "4qy1_324", "4qy2_324" ],

"_A_ H16N3 A/B/C" => [ "4f23_329", "4fiu_329" ],

"_A_ H17N10 A/B C/D E/F" => [ "4h32_335" ],

"_A_ H18N11 A/B/C" => [ "4mc5_330" ],

"_B_ Victoria A/B C/D E/F" => [ "4fqm_348" ]
);



my $workdir = "$path_environment/Interactions/Data/PDB_splitable/Trimeric";
my $outdir = "$path_environment/Interactions/Data/PDB_Clean/Trimeric";


#Foreach chain name combinations
foreach $_ (keys %HashPDBinfo) {


  my $Flu; #A ou B
  my $State; #Rien = Prefusion, sinon Intermediate or PostFusion
  my $Subtype;

  my $Namechain1; #HA1 Proto1
  my $Namechain2; #HA2 Proto1
  my $Namechain3; #HA1 Proto2
  my $Namechain4; #HA2 Proto2
  my $Namechain5; #HA1 Proto3
  my $Namechain6; #HA2 Proto3



  my @SplitScalarChain = split (" ", $_);
    #EXEMPLES of $_ :
    # CASE 1 : "_A_ H1N1 A/D/F";
    # CASE 3 : "_A_ H1N1 A/D/F"; with Postfusion or Intermediates between subtype and A/B

    # CASE 2 : "_A_ H3N2 A/B C/D E/F";
    # CASE 4 : "_A_ H3N2 A/B C/D E/F"; with Postfusion or Intermediates between subtype and A/B


  # IF CASE 1
  if (scalar (@SplitScalarChain) == 3) {

      $_ =~ m/_([AB])_\s(.+)\s(.)\/(.)\/(.)/;

      $State = "Prefusion"; $Flu = $1; $Subtype = $2;
      $Namechain1 = $3; $Namechain2 = $3;   $Namechain3 = $4; $Namechain4 = $4; $Namechain5 = $5; $Namechain6 = $5;
      print "$1 $2 $3/$4/$5\n";
      }


  # IF CASE 2
  elsif (scalar (@SplitScalarChain) == 4) {

      $_ =~ m/_([AB])_\s(.+)\s(.+)\s(.)\/(.)\/(.)/;

      $Flu = $1; $State = $2; $Subtype = $3;
      $Namechain1 = $4; $Namechain2 = $4;
      $Namechain3 = $5; $Namechain4 = $5;
      $Namechain5 = $6; $Namechain6 = $6;
      print "$1 $2 $3 $4/$5/$6\n";
    }


   # IF CASE 3
   elsif (scalar (@SplitScalarChain) == 5) {

      $_ =~ m/_([AB])_\s(.+)\s(.)\/(.)\s(.)\/(.)\s(.)\/(.)/;

      $State = "Prefusion"; $Flu = $1; $Subtype = $2;
      $Namechain1 = $3; $Namechain2 = $4;
      $Namechain3 = $5; $Namechain4 = $6;
      $Namechain5 = $7; $Namechain6 = $8;
      print "$1 $2 $3/$4 $5/$6 $7/$8\n";
      }


   # IF CASE 4
   elsif (scalar (@SplitScalarChain) == 6) {

      $_ =~ m/_([AB])_\s(.+)\s(.+)\s(.)\/(.)\s(.)\/(.)\s(.)\/(.)/;

      $Flu = $1; $State = $2; $Subtype = $3;
      $Namechain1 = $4; $Namechain2 = $5;
      $Namechain3 = $6; $Namechain4 = $7;
      $Namechain5 = $8; $Namechain6 = $9;
      print "$1 $2 $3 $4/$5 $6/$7 $8/$9\n";
      }


    #INTITIALIZE OUTPUT PATH AND FILE
    my $dir_output = "$outdir/Flu$Flu/$State/$Subtype";
    make_path($dir_output) unless -d $dir_output;



    #Foreach PDB file associated to the KEY
    foreach $_ (@{$HashPDBinfo{$_}}) {


      #GET THE NUMBER OF FIRST HA2 RESIDUE
      $_ =~ m/(.+)_([0-9]+)/;
      my $CodePDB = $1;   my $NumberFirstResidueHA2 = $2;


      open (FICHIERPDB	  , "<$workdir/Flu$Flu/$State/$Subtype/$CodePDB.pdb") or die "Pb d'ouverture : $! $Flu/$State/$Subtype/$CodePDB";
      open (FICHIERPDBCLEAN ,">$dir_output/PDBClean$CodePDB.txt") or die "Pb d'ouverture : $! $dir_output/$CodePDB";




      my (@HA1_PROTO1, @HA2_PROTO1, @HA1_PROTO2, @HA2_PROTO2, @HA1_PROTO3, @HA2_PROTO3);
      my @ArrayOfArrays = (\@HA1_PROTO1, \@HA2_PROTO1, \@HA1_PROTO2, \@HA2_PROTO2, \@HA1_PROTO3, \@HA2_PROTO3);

     #LOOP OVER PDB FILE
      while ( my $line = <FICHIERPDB> )
      {
        if ( $line =~ m/^ATOM/)
        {

	         my @list=split(" ",$line);
          #ATOM    138  O   VAL A  23      35.057  66.370 122.942  1.00 33.18           O

           #UNLESS ATOM = HYDROGEN
           unless ($list[2] =~ m/^H.*/) {


          #CLEAN/HOMOGENEIZE PDB FILES WITH NEW CHAIN NAME NOMENCLATURE + NEW NUMERATION


           #IF NUMBER RESIDUE < Number first residu HA2 --> WE ARE IN HA1
       		if ($list[4] eq "$Namechain1" and $list[5] < $NumberFirstResidueHA2) {
            my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto1", $list[5]);
            push (@HA1_PROTO1, $JoinInfos);

       	}

          # IF NUMBER RESIDUE >= Number first residu HA2 --> WE ARE IN HA2
          # We change the previous number of the residu by the new number (Numeration HA2 starting from 1)
       		elsif ($list[4] eq "$Namechain2" and $list[5] >= $NumberFirstResidueHA2) {

             #NEW HA2 NUMERATION
       		   my $NumberClean = $list[5] - $NumberFirstResidueHA2 + 1;

             my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $NumberClean, $list[6], $list[7], $list[8], "HA2", "Proto1", $list[5]);
             push (@HA2_PROTO1, $JoinInfos);

       	}




       		elsif ($list[4] eq "$Namechain3" and $list[5] < $NumberFirstResidueHA2) {
            my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto2", $list[5]);
            push (@HA1_PROTO2, $JoinInfos);

       	}


       		elsif ($list[4] eq "$Namechain4" and $list[5] >= $NumberFirstResidueHA2) {

             #NEW HA2 NUMERATION
       		   my $NumberClean = $list[5] - $NumberFirstResidueHA2 + 1;

             my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $NumberClean, $list[6], $list[7], $list[8], "HA2", "Proto2", $list[5]);
             push (@HA2_PROTO2, $JoinInfos);

      }





       		elsif ($list[4] eq "$Namechain5" and $list[5] < $NumberFirstResidueHA2) {
            my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto3", $list[5]);
            push (@HA1_PROTO3, $JoinInfos);

       	}

       		elsif ($list[4] eq "$Namechain6" and $list[5] >= $NumberFirstResidueHA2) {

             #NEW HA2 NUMERATION
       		   my $NumberClean = $list[5] - $NumberFirstResidueHA2 + 1;

             my $JoinInfos = join " ", ($list[1], $list[2], $list[3], $list[4], $NumberClean, $list[6], $list[7], $list[8], "HA2", "Proto3", $list[5]);
             push (@HA2_PROTO3, $JoinInfos);

       	}

            }


    }
}


#PRINT FILE
my $DoColumns = "%-10s\t"x10 . "%-10s\n";
foreach my $Array (@ArrayOfArrays) {
  foreach my $InfosLine (@$Array) {

  my @SplitForColumns = split (" ", $InfosLine);
  printf FICHIERPDBCLEAN $DoColumns, @SplitForColumns;

  }
}



  close FICHIERPDB;
  close FICHIERPDBCLEAN;
      }

		       }
