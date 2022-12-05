use warnings;
use strict;

use File::Path qw( make_path );
use List::Util qw(sum uniq any max reduce);
use JSON;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;


use List::Compare;



#CREATED FASTA FORMAT FOREACH CLEANED FILE IN ORDER TO FEED 3DCOFFEE (Sequence aligment based on structural informations)
#WE CREATE A FASTA FORMAT PER PROTOMER OF A SIGNLE STRUCTURE (Because sometimes some amino acids are not resolved for a single protomer)


my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



my $workdir = "$path_environment/Interactions/Data/PDB_Clean/Trimeric";
my $outdir = "$path_environment/Interactions/Data/FastaFiles";
my $out_tcoffee = "$path_environment/Interactions/Data/Templates_tcoffee";




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


my $Nb = 0;
my $NbTemplateTakenAndRejected = 0;
my $JSONdata;
my @CheckCodePDBProcedes = ();



open(ARRAYFILE, "<", "$path_environment/Interactions/Data/ArrayCleanWithRefTableauPDB.dat");
$JSONdata = <ARRAYFILE>; close(ARRAYFILE);
my $ArrayRefTableauPDBClean = decode_json($JSONdata);


my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");
my @Protos = ("Proto1", "Proto2", "Proto3");


foreach my $flu (@Flu) {

  open (OUTPUTFASTAALL, ">$outdir/$flu/FastaFileAll_$flu.fasta") or die "Unable to open FastaFileAllPDB.fasta:$!\n";

  #FASTA template file that will be used for sequence alignment using structural informations (3Dcofee)
  open (TEMPLATEFORTCOFFEE, ">$out_tcoffee/Template3DCoffee_$flu.txt") or die "Unable to open Template3DCoffee_$flu.txt:$!\n";
  my $NbTemplateTcoffee = 0;
  my $NbrejectedStructuresForTemplateCoffee = 0;



  foreach my $StructureState (@StructureStates){

    unless ($flu eq "FluB" and $StructureState eq "Intermediates" or $flu eq "FluB" and $StructureState eq "Postfusion"){


    opendir(DIR, "$workdir/$flu/$StructureState") or die "Unable to enter dir $workdir/$flu/$StructureState:$!\n";
    my @FolderNamesSubtypes = readdir(DIR) or die "Unable to read $workdir/$flu/$StructureState:$!\n";


    LOOP:foreach my $Subtype (@FolderNamesSubtypes){

        next if ($Subtype =~ m/^\./);
        next if ($Subtype =~ m/^\./);

        if (-d "$workdir/$flu/$StructureState/$Subtype"){                  # is this a directory?

        my $dir_output = "$outdir/$flu/$StructureState/$Subtype";
        make_path($dir_output);
        open (OUTPUTFASTA, ">$dir_output/FastaFile$Subtype.fasta") or die "Unable to open FastaFile$Subtype.fasta:$!\n";

        opendir(DIRTWO, "$workdir/$flu/$StructureState/$Subtype") or die "Unable to open $Subtype:$!\n";
        my @filesCodePDB = readdir(DIRTWO);



            foreach my $CodePDB (@filesCodePDB){


              next if ($CodePDB =~ m/^\./);
              next if ($CodePDB =~ m/^\./);

              if ($CodePDB =~ m/^PDBClean(.+)\./) {push @CheckCodePDBProcedes, $1; $Nb++; }

                unless (-d "$workdir/$flu/$StructureState/$Subtype/$CodePDB") {


                  open (CODEPDB, "$workdir/$flu/$StructureState/$Subtype/$CodePDB") or die "Unable to open $CodePDB:$!\n";

                  #INITIALISE NEW VARIABLE FOR EACH PDB CODE
                  my @Grep = ();
                  my $InfosCodePDB;
                  my $UpperCaseCODPDB;
                  my $LowerCaseCODPDB;
                  my @Keys;

                  #Grep Code PDB IN array containing all PDB codes and associated informations
                  if ($CodePDB =~ m/PDBClean(.+)\.txt/) { $UpperCaseCODPDB = uc($1);
                                                        $LowerCaseCODPDB = lc($1);
                                                        @Grep = grep { m/^$UpperCaseCODPDB/ } @$ArrayRefTableauPDBClean;

                       }


                       ##### WARNINGS CALLS IF SOMETHING UNEXPECTED
                       my $ScalarGrep = scalar (@Grep);

                       print Dumper \@Grep if $ScalarGrep == 0;
                       print "$CodePDB   NO GREP NO CAN'T PROCEED\n\n" if $ScalarGrep == 0;

                       if (scalar (@Grep) >= 2) { print "WARNING GREP $CodePDB MORE THAN ONE FOUND\n\n"; last; }



                       foreach (@Protos) {

                         $InfosCodePDB = $Grep[0]; $InfosCodePDB =~ m/^($UpperCaseCODPDB\|.{4,6})\|.+/;
                         $InfosCodePDB = $LowerCaseCODPDB . $_;

                         ## NECESARRY PROCESS TO MAKE 3DCOFFEE TEMPLATE FILE
                         my $PathPDB_ForTemplate = "$path_environment/Interactions/Data/PDB/Trimeric/$flu/$StructureState/$Subtype/$LowerCaseCODPDB.pdb";
                         $NbTemplateTcoffee++, print TEMPLATEFORTCOFFEE "\n>$InfosCodePDB _P_ $PathPDB_ForTemplate\n" unless $StructureState ne "Prefusion";
                         $NbrejectedStructuresForTemplateCoffee++ if $StructureState =~ /Postfusion|Intermediates/;

                       }

                       #Make Fasta files thanks to opened cleaned pdb file
                       CleanToFastaAllProto (\%ThreeLetterToOneLetterAACode, $LowerCaseCODPDB);

              close (CODEPDB);
                }
                                  }
                      close OUTPUTFASTA;   }

             else {print "$workdir/$flu/$StructureState/$Subtype not a directory so can't run\n"; next LOOP; }

    close (DIRTWO);
#        close NBPERPROTO;
      }
closedir(DIR);
}
 }
close OUTPUTFASTAALL; close TEMPLATEFORTCOFFEE;


$NbTemplateTcoffee=$NbTemplateTcoffee/3;
$NbrejectedStructuresForTemplateCoffee=$NbrejectedStructuresForTemplateCoffee/3;

$NbTemplateTakenAndRejected += $NbTemplateTcoffee;
$NbTemplateTakenAndRejected += $NbrejectedStructuresForTemplateCoffee;
## DIVIDE BY 3 CAR CALCUL FOR 3 PROTO...


print "Template alignement T-coffee = $NbTemplateTcoffee for $flu\nStructures not considered cuz postfusion or intermediates $NbrejectedStructuresForTemplateCoffee for $flu\n\n";
}

print "Structure template taken + not considered = $NbTemplateTakenAndRejected\n\n";





# CHECK IF ALL FILES HAVE BEEN PROCESSED
my @PDBArrayClean;
foreach $_ (@$ArrayRefTableauPDBClean) { my @Split = split ('\|', $_); push @PDBArrayClean, $Split[0];  }

@CheckCodePDBProcedes = map { uc } @CheckCodePDBProcedes;

my $CompareArray = List::Compare->new( \@CheckCodePDBProcedes, \@PDBArrayClean );
my @WeDidntMakeThisFile = $CompareArray->get_Lonly;
unless (scalar (@WeDidntMakeThisFile) == 0) {
  foreach $_ (@WeDidntMakeThisFile) {print " WARNING WE DON'T HAVE THIS FILE IN SECOND ARRAY : $_\n";}
  }

#print Dumper \@CheckCodePDBProcedes;
#print "\n\n\n\n\n";
#print Dumper \@PDBArrayClean;
#print "\n\n\n\n\n";

print "Sequences processed = $Nb\n";







my $NumberSeqFastaAll = 0;

foreach my $flu (@Flu) {
open (FASTAALL, "<$outdir/$flu/FastaFileAll_$flu.fasta") or die "Unable to open FastaFileAllPDB.fasta:$!\n";

while ( $_ = <FASTAALL> )
{
$NumberSeqFastaAll++ if $_ =~ m/^>.+/;
  }
    }
    $NumberSeqFastaAll=$NumberSeqFastaAll/3; ##Cuz calcul for 3 protos
print "Fasta sequences created (fluA + fluB) = $NumberSeqFastaAll\n";







####################################################################






sub CleanToFastaAllProto {
  my ($ref_HashConvertion, $LowerCaseCODPDB) = @_;


  my @Protomers = ("Proto1", "Proto2", "Proto3");
  my @Lines_pdb= ();



  my @Previous = ();
  #Put file in an array
  while (my $line = <CODEPDB>){
    push @Lines_pdb, $line;
      }



    foreach my $Protomer (@Protomers) {


      print OUTPUTFASTA "\n"; print OUTPUTFASTAALL "\n";


      my $Infos = $LowerCaseCODPDB . $Protomer;

      #Print Fasta header
      print OUTPUTFASTA ">$Infos\n"; print OUTPUTFASTAALL ">$Infos\n";

      LOOP:foreach my $LineCodePDB (@Lines_pdb) {
          my @SplitLine = split (" ", $LineCodePDB);
          #  5 CB ASP A 5 34.370 88.449 159.751 HA1 Proto1


          if (($SplitLine[8] eq "HA1" and $SplitLine[9] eq "$Protomer") or ($SplitLine[8] eq "HA2" and $SplitLine[9] eq "$Protomer")) {

            # IF / ELSIF / ELSO ENABLING TO PRINT ONE RESIDU ONLY NO MATTER HOW MANY ATOMS ARE DESCRIBED
            if( scalar(@Previous) == 0) { print OUTPUTFASTA "$ref_HashConvertion->{$SplitLine[2]}"; print OUTPUTFASTAALL "$ref_HashConvertion->{$SplitLine[2]}";}
            elsif ( $Previous[0] eq $SplitLine[2] and $Previous[1] eq $SplitLine[4]) { next LOOP;}
            else {

              print OUTPUTFASTA "$ref_HashConvertion->{$SplitLine[2]}"; print OUTPUTFASTAALL "$ref_HashConvertion->{$SplitLine[2]}";

                }
              }

  @Previous = ($SplitLine[2], $SplitLine[4]);

#print "$Infos\n";


   }
 }
}
