use warnings;
use strict;

use File::Path qw( make_path );




my $path_environment = "/Users/valentin.ozl/GitHub";
#    Change your environment : replace "/Users/valentin.ozl/GitHub" by your work environment



my @Flu = ("FluA", "FluB");
my @StructureStates = ("Prefusion", "Intermediates", "Postfusion");
my $NbTotalStructures = 0;


#Foreach influenza virus type
foreach my $flu (@Flu) {
  my $StructureFlu = 0;



  #Foreach state of structure : prefusion, intermediary, postufusion
  foreach my $StructureState (@StructureStates){
    my $StructureFluState = 0;


    #Initiate work and output DIR
    my $workdir = "$path_environment/Interactions/Data/PDB/Trimeric/$flu/$StructureState";
    my $outdir = "$path_environment/Interactions/Data/PDB_splitable/Trimeric/$flu/$StructureState";

    opendir(DIR, "$workdir") or die "Unable to enter dir $workdir:$!\n";
    my @FolderNamesSubtypes = readdir(DIR) or die "Unable to read $workdir:$!\n";

    #Grab all subtypes
    @FolderNamesSubtypes =  grep { ! m/^\./ } @FolderNamesSubtypes;
    my $NbSubtypes = scalar (@FolderNamesSubtypes);



    #Foreach grabbed subtype
     foreach my $Subtype (@FolderNamesSubtypes) {

       my $NbStructuresPerSubtype;
        next if ($Subtype eq ".");
        next if ($Subtype eq "..");


        # is this an existing directory?
        if (-d "$workdir/$Subtype"){

        opendir(DIRTWO, "$workdir/$Subtype") or die "Unable to open $workdir/$Subtype:$!\n";
        my @filesCodePDB = readdir(DIRTWO);

        #Grab number structure per subtype
        @filesCodePDB = grep { ! m/^\./ } @filesCodePDB; #On ne comptabilise pas les fichiers cachés
        my $NbStructuresPerSubtype = scalar (@filesCodePDB);

        print "$flu $StructureState // Subtypes $Subtype, number structures processed : $NbStructuresPerSubtype\n";

        my $dir_outputSubtype = "$outdir/$Subtype";
        make_path($dir_outputSubtype);



        #Foreach structure
        foreach my $CodePDB (@filesCodePDB){
          $NbTotalStructures++;
          $StructureFlu++;
          $StructureFluState++;

          next if ($CodePDB eq ".");
          next if ($CodePDB eq "..");

          #Open PDB file input
          open (READCODEPDB, "<$workdir/$Subtype/$CodePDB") or die "Unable to open $workdir/$Subtype/$CodePDB:$!\n";

          my @file_lines = <READCODEPDB>;
          close(READCODEPDB);

          open( WRITEFILE, ">$dir_outputSubtype/$CodePDB" ) or die "Unable to open $workdir/$Subtype/$CodePDB:$!\n";



          INNERFOR:for (my $i = 0; $i <= $#file_lines; $i++) {


            #COLUMNS DATA TYPE FIELD DEFINITION
            #-------------------------------------------------------------------------------------
            #ATOM    342  CG  LEU A  54      28.908 -19.234  47.872  1.00 28.03           C

            #1 - 6 Record name "ATOM "
            #7 - 11 Integer serial Atom serial number.
            #13 - 16 Atom name Atom name.
            #17 Character altLoc Alternate location indicator.
            #18 - 20 Residue name resName Residue name.
            #22 Character chainID Chain identifier.
            #23 - 26 Integer resSeq Residue sequence number.
            #27 AChar iCode Code for insertion of residues.
            #31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms.
            #39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms.
            #47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms.
            #55 - 60 Real(6.2) occupancy Occupancy.
            #61 - 66 Real(6.2) tempFactor Temperature factor.
            #77 - 78 LString(2) element Element symbol, right-justified.
            #79 - 80 LString(2) charge Charge on the atom.


            #Read line containing ATOM
            if ($file_lines[$i] =~ /ATOM  (.{5}).{1}(.{4})(.{1})(.{3}).{1}(.{1})(.{4})(.{1}).{3}(.{8})(.{8})(.{8})(.{6})(.{6}).{10}(.{2})(.{2})/)
            {

              #$AtomNb = $1; $AtomName = $2;$ResAltLoc = $3;$AAname = $4; $ChainId = $5;
              #$AAnb = $6; $insertion = $7; $x = $8; $y = $9; $z = $10; $occupancy = $11;
              #$TempFactor = $12; $elementSymbol = $13; $atomCharge = $14;

              print WRITEFILE "ATOM $1 $2 $4 $5 $6 $8 $9 $10 $11 $12 $13 $14\n";




            }



                                                                    }
               close (WRITEFILE);
       }
     }
   }
print "$flu $StructureState --> Total Subtypes processed : $NbSubtypes // Total Structures processed : $StructureFluState\n\n";
  }
print "$flu $StructureFlu structures processed\n\n\n";
}
print "Total structures processed : $NbTotalStructures\n\n";
