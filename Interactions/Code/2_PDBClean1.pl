use warnings;
use strict;

use File::Path qw( make_path );

############# CODE PDB CLEAN for structures which delimit HA1 and HA2 but has non-conventionnal chain names to identify HA1 and HA2 of a certain protomer
############# ALSO HOMOGENEIZE AA NUMERATION (Some files have HA1 : 1-325 and HA2 1-150 (approximatively) BUT SOME have HA1 : 1-325 and HA2 : 325-500)


my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



#######################

# Hash of all identified chain names combinations (KEYS) within raw PDB files and corresponding PDB files (Value = ref list)

my %HashPDBinfo = (

"_A_ H1N1 A/B C/D E/F" => [ "1ru7", "1rvx", "1rvz", "3al4", "3lzg", "3m6s", "3ube", "3ubj", "3ubn",
"3ubq", "4edb", "4eef", "4gxu", "4gxx", "4jtv", "4jtx", "4ju0", "4jug", "4juh", "4juj",
"4lxv", "4m4y", "5ugy", "5ujz", "5uk0", "5uk1", "5uk2", "5vmc", "5vmf", "5vmg", "5vmj",
"6hjn", "6hjp", "6hjq", "6hjr", "1rd8", "6wj1", "6xgc" ],

"_A_ H1N1 H/I J/K L/M" => [ "1ruy", "1ruz", "2wrg", "2wrh", "1rv0", "1rvt" ],

"_A_ H1N1 C/D E/F G/H" => [ "5w6r" ],

"_A_ H1N1 M/P N/Q O/R" => [ "5wko" ],

"_A_ H1N1 A/B C/D I/J" => [ "6lks" ],

"_A_ H1N2 A/B C/D E/F" => [ "4f3z" ],

"_A_ H2N2 A/B C/D E/F" => [ "3qqo", "4hlz" ],

"_A_ H2N2 A/a B/b C/c" => [ "4hg4" ],

"_A_ H2N3 A/B C/D E/F" => [ "4w8n" ],

"_A_ H3N2 A/B C/D E/F" => [ "1hgd", "1hge", "1hgf", "1hgg", "1hgh", "1hgi", "1hgj", "1ken", "2hmg",
"2ypg", "3hmg", "3vun", "3ztj", "4fnk", "4fqr", "4gms", "4hmg", "4nm8", "4o5i", "4zcj",
"5hmg", "5k9q", "5kan", "5t6n", "5thf", "5vtq", "5vtr", "5vtu", "5vtv", "5vtw",
"5vtx", "5vty", "5vtz", "5vu4", "5xrs", "5xrt", "6bkm", "6bkq", "6cex", "6tzb", "3eym",
"6nhp", "6nhq", "6nhr"  ],

"_A_ Intermediates H3N2 A/B C/D E/F" => [ "6y5g", "6y5h", "6y5i", "6y5j", "6y5k", "6y5l", "7k39", "7k37", "7k3a", "7k3b"  ],

"_A_ H3N2 A/B C/D G/H" => [ "6n4f" ],

"_A_ H3N2 A/B E/F I/J" => [ "4ubd" ],

"_A_ H3N8 A/B C/D E/F" => [ "4uo0", "4uo1", "4uny", "4unz", "4unw", "4uo3", "4uo7", "4uo5", "4unx" ],

"_A_ H3N8 A/B D/E G/H" => [ "1mql", "1mqn", "1mqm" ],

"_A_ H4N6 A/B C/D E/F" => [ "6v44" ],

"_A_ H5N1 A/B C/D E/F" => [ "2fk0", "2ibx", "3fku", "3s11", "4bgz", "4bh0", "4bh1", "4jul", "4jun",
"4kdq", "4kth", "4mhh", "4mhi", "5e2y", "5e2z", "5e30", "5jw4", "6cf5", "6e7g", "6e7h",
"3znk", "3znl", "3znm", "4cqw", "4cqx", "4cqy", "4cqv", "6vmz" ],

"_A_ H5N1 A/B C/D G/I" => [ "4mhj" ],

"_A_ H5N1 A/B G/H I/J" => [ "4n5y", "4n5z" ],

"_A_ H5N1 A/B G/L K/H" => [ "6b3m" ],

"_A_ H5N6 A/B C/D E/F" => [ "5hu8" ],

"_A_ H5N8 A/B C/D E/F" => [ "5huf" ],

"_A_ H6N1 A/B C/D E/F" => [ "4xkd", "4xke", "4xkf", "4xkg", "5t0b", "5t0d", "5t0e", "5t08",
 "4wsv", "4wsu" ],

"_A_ H6N1 A/B G/H I/J" => [ "4wst" ],

"_A_ H6N6 A/B C/D E/F" => [ "5bny", "5bqy", "5bqz" ],

"_A_ H7N2 A/B C/D E/F" => [ "3m5g", "3m5h", "3m5i", "3m5j" ],

"_A_ H7N2 A/C B/D H/J" => [ "6mlm" ],

"_A_ H7N7 A/B C/D E/F" => [ "4dj6", "4dj7", "4dj8", "4fqv" ],

"_A_ H7N9 A/B C/D E/F" => [ "4ln3", "4ln4", "4ln6", "4ln8", "5t6s", "6d7u", "6d8b", "6d8d", "6ii9" ],

"_A_ H7N9 A/B E/F I/J" => [ "6d7c" ],

"_A_ H7N9 A/B D/E G/H" => [ "6fyu" ],

"_A_ H8N4 A/B C/D E/F" => [ "6v46" ],

"_A_ H10N2 A/B C/D E/F" => [ "4cyv", "4cyw", "4cyz", "4cz0" ],

"_A_ H10N7 A/B C/D E/F" => [ "4wsw", "6tvt", "6tjw", "6txo", "6tvc", "6tva", "6tvb", "6twh", "6twi", "4d00" ],

"_A_ H10N7 C/D E/F K/L" => [ "6tvs", "6twv" ],

"_A_ H10N7 A/B G/H I/J" => [ "6tvr", "6tjy", "6tvd", "6tvf" ],

"_A_ H10N7 A/B E/F I/J" => [ "6ty1" ],

"_A_ H10N7 A/B G/H C/D" => [ "6tws" ],

"_A_ H10N8 A/B C/D E/F" => [ "4xq5", "4xqo", "4xqu", "5tgo", "5tgu", "5tgv", "5th0", "5th1", "5thb",
 "5thc" ],

"_A_ H10N8 A/B O/P W/X" => [ "4wsx" ],

"_A_ H11N9 A/B C/D E/F" => [ "6v47" ],

"_A_ H13N6 A/B C/D E/F" => [ "4kpq", "4kps" ],

"_A_ H14N5 A/B E/F K/L" => [ "6v48" ],

"_A_ H15N9 A/B C/D E/F" => [ "6v49" ],

"_A_ H18N11 A/B C/D E/F" => [ "4k3x" ],

"_A_ Postfusion H3N2 -/B -/D -/F" => [ "1htm" ],

"_A_ Postfusion H3N2 -/A -/B -/C" => [ "1qu1" ],

"_B_ Ancestral A/B C/D E/F" => [ "4nrj", "4nrk", "4nrl" ],

"_B_ Yamagata A/B C/D E/F" => [ "4m40", "4m44"]

);





my $workdir = "$path_environment/Interactions/Data/PDB_splitable/Trimeric";
my $outdir = "$path_environment/Interactions/Data/PDB_Clean/Trimeric";




#Foreach chain name combinations
foreach $_ (keys %HashPDBinfo) {

  my $Flu; #A ou B
  my $State; #Nothing = Prefusion, otherwise Intermediate or PostFusion
  my $Subtype;

  my $Namechain1; #HA1 Proto1
  my $Namechain2; #HA2 Proto1
  my $Namechain3; #HA1 Proto2
  my $Namechain4; #HA2 Proto2
  my $Namechain5; #HA1 Proto3
  my $Namechain6; #HA2 Proto3


  my @SplitScalarChain = split (" ", $_);
  #EXEMPLES of $_ :
  # CASE 1 : "_A_ H1N1 A/B C/D E/F"
  # CASE 2 : "_A_ Postfusion H1N1 A/B C/D E/F"

  #REDEFINE THE CHAIN NAMES TO HOMOGENEIZE HA1/HA2 NOMENCLATURE

  # If case 1:
  if (scalar (@SplitScalarChain) == 5) {

      $_ =~ m/_([AB])_\s(.+)\s(.)\/(.)\s(.)\/(.)\s(.)\/(.)/;

      #Identify and asign new variable
      $State = "Prefusion"; $Flu = $1; $Subtype = $2;
      $Namechain1 = $3; $Namechain2 = $4;
      $Namechain3 = $5; $Namechain4 = $6;
      $Namechain5 = $7; $Namechain6 = $8;

      print "$1 $2 $3 $4 $5 $6 $7 $8\n";
      }

  # If case 2:
  elsif (scalar (@SplitScalarChain) == 6) {

      $_ =~ m/_([AB])_\s(.+)\s(.+)\s(.)\/(.)\s(.)\/(.)\s(.)\/(.)/;

      #Identify and asign new variable
      $Flu = $1; $State = $2; $Subtype = $3;
      $Namechain1 = $4; $Namechain2 = $5;
      $Namechain3 = $6; $Namechain4 = $7;
      $Namechain5 = $8; $Namechain6 = $9;
      print "$1 $2 $3 $4 $5 $6 $7 $8 $9\n";
    }



  #Initialize OUTPUT path and file
  my $dir_output = "$outdir/Flu$Flu/$State/$Subtype";
  make_path($dir_output) unless -d $dir_output;



  #Foreach PDB file belonging to the KEY
  foreach my $CodePDB (@{$HashPDBinfo{$_}}) {

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



              #CLEAN/HOMOGENEIZE PDB FILES WITH NEW NOMENCLATURE
              #PROTOMER 1
	            if ($list[4] eq $Namechain1) {
                  my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto1", $list[5]);
                  push (@HA1_PROTO1, $JoinInfos);

				            }
	            elsif ($list[4] eq $Namechain2 ) {
			              if ($list[5] >= 500) {

                        #ENABLE NUMERATION STARTING FROM 1 WHEN WE ARE IN HA2
		                    my $CleanNumber = $list[5] - 500;

                        my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $CleanNumber, $list[6], $list[7], $list[8], "HA2", "Proto1", $list[5]);
                        push (@HA2_PROTO1, $JoinInfos);

			                }
  		              else {
                      my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA2", "Proto1", $list[5]);
                      push (@HA2_PROTO1, $JoinInfos);

					                 }
				            }





              #PROTOMER 2
	            elsif ($list[4] eq $Namechain3 ) {
                my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto2", $list[5]);
                push (@HA1_PROTO2, $JoinInfos);

	               }

	            elsif ($list[4] eq $Namechain4 ) {
					        if ($list[5] >= 500) {

                     #ENABLE NUMERATION STARTING FROM 1 WHEN WE ARE IN HA2
							       my $CleanNumber = $list[5] - 500;
                     my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $CleanNumber, $list[6], $list[7], $list[8], "HA2", "Proto2", $list[5]);
                     push (@HA2_PROTO2, $JoinInfos);

							      	}

					        else {
                     my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA2", "Proto2", $list[5]);
                     push (@HA2_PROTO2, $JoinInfos);

					              }
				              }





              #PROTOMER 3
	            elsif ($list[4] eq $Namechain5 )	{
                my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA1", "Proto3", $list[5]);
                push (@HA1_PROTO3, $JoinInfos);

	              }

	            elsif ($list[4] eq $Namechain6 ) {

                  if ($list[5] >= 500) {

                     #ENABLE NUMERATION STARTING FROM 1 WHEN WE ARE IN HA2
							       my $CleanNumber = $list[5] - 500;
                     my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $CleanNumber, $list[6], $list[7], $list[8], "HA2", "Proto3", $list[5]);
                     push (@HA2_PROTO3, $JoinInfos);

								      }

					        else {
                     my $JoinInfos = join "\t", ($list[1], $list[2], $list[3], $list[4], $list[5], $list[6], $list[7], $list[8], "HA2", "Proto3", $list[5]);
                     push (@HA2_PROTO3, $JoinInfos);

					               }
				             }
          }


  }
}

  #PRINT CLEANED FILE
  my $DoColumns = "%-10s\t"x10 . "%-10s\n";
  foreach my $Array (@ArrayOfArrays) {
    foreach my $InfosLine (@$Array) {

    my @SplitForColumns = split ("\t", $InfosLine);
    printf FICHIERPDBCLEAN $DoColumns, @SplitForColumns;

    }
}



  close FICHIERPDB;
  close FICHIERPDBCLEAN;
      }

		       }
