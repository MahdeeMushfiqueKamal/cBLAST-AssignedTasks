#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use File::Find;
use File::Path;
use Cwd;
#use diagnostics;

 ########################### WELCOME  ##############################
#                                                                   #
#                          sgRNAcas9                                # 
# ---a tool for fast designing CRISPR sgRNA with high specificity   #
#                                                                   #
# AUTHOR  : Xie Shengsong                                           #
# Email   : ssxieinfo\@gmail.com                                    #
# Homepage: www.biootools.com                                       #
#           BiooTools (Biological online tools)                     #  
# Shanghai Institute of Biochemistry and Cell Biology (SIBCB)       #
# School of Veterinary Medicine, Huazhong Agricultural University   #
# Version: sgRNAcas9_3.0.5                                          #
# Begin       : 2013.12.9                                           #
# LAST REVISED: 2015.7.28                                           #
 ###################################################################

my $oldtime = time();

my ($Inputfile_Fasta, $truncat, $GC_l, $GC_m, $Genome, $Option, $Type, $Seqmap_vesion, $Num_mismatch, $offset_s, $offset_e, $path);

GetOptions( "i=s" => \$Inputfile_Fasta,        #Input file
            "x=i" => \$truncat,                #Length of sgRNA[20]
            "l=i" => \$GC_l,                   #The minimum value of GC content [20]
			"m=i" => \$GC_m,                   #The maximum value of GC content [80] 
	        "g=s" => \$Genome,                 #The reference genome sequence
			"o=s" => \$Option,                 #Searching CRISPR target sites using DNA strands based option(s/a/b)
			"t=s" => \$Type,                   #Type of gRNA searching mode(s/p) 
			"v=s" => \$Seqmap_vesion,          #Operation system [w, for windows; l, for linux-64; u, for linux-32;  m, for MacOSX-64; a, for MacOSX-32]
			"n=i" => \$Num_mismatch,           #Maximum number of mismatches [5]
			"s=i" => \$offset_s,               #The minimum value of sgRNA offset [-2]
            "e=i" => \$offset_e,               #The maximum value of sgRNA offset [32]
            "p=s" => \$path,                   #Output path
          );

#default
$truncat ||= "20";  
$GC_l ||= "20";                                #20 % < GC% < 80 %
$GC_m ||= "80";
$Seqmap_vesion ||= "l";                        #linux
$Num_mismatch ||="5";                          #number of mismatches: 5
$offset_s ||="-3";                             #sgRNA offset: -2 to 32 bp
$offset_e ||="33";
my $dir_default = getcwd;                      #default output
$path ||= $dir_default;

my $dir =$path;

mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta",0755)||die "Can't create directory: Directory exists at $dir. Please delete, move or rename the exist directory before you run this program.$!" ;

open  (LOG, ">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/sgRNAcas9.Log.txt") || die "Can't open sgRNAcas9.Log.txt for writing!" ."\n";
print  LOG "################################# Log ###########################################".        "\n\n";
#print "Writing Log information.                                                                      " ."\n";
print  LOG "#                              sgRNAcas9                                                  " ."\n";
print  LOG "#     ---a tool for fast designing CRISPR sgRNA with high specificity                     " ."\n";          
print  LOG "#                                                                                         " ."\n";
print  LOG "#       contact:  Xie Shengsong, Email: ssxieinfo\@gmail.com                                .\n\n";
######


#mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/A.Final_report",0755)||die "can't create directory: $!" ;

print "\n\tWelcome to sgRNAcas9\n";
print "\t---a tool for fast designing CRISPR sgRNA with high specificity\n";
print "\t---------------------------------------------------------\n";
print "Version   : 3.0"."\n";
print "Copyright : Free software"."\n";
print "Author    : Shengsong Xie"."\n";
print "Email     : ssxieinfo\@gmail.com"."\n";
print "Homepage  : www.biootools.com"."\n";

my $local_time;
$local_time = localtime();
print "Today     : $local_time\n\n";
print  LOG "# Time, begin at $local_time."."\n";
print  LOG "# Usage: perl $0 -i $Inputfile_Fasta -x $truncat -l $GC_l -m $GC_m -g $Genome -o $Option -t $Type -v $Seqmap_vesion -n $Num_mismatch -s $offset_s -e $offset_e -p $path 2>log.txt\n\n";

################################### format seq ###################################
print "Start sgRNAcas9 program........\n";
print "Step1: Format target sequences.\n";
print  LOG "# Start sgRNAcas9 program........\n";
print  LOG "# Step1: Format target sequences.\n";

open (Inseq, $Inputfile_Fasta) || die "Can't open $Inputfile_Fasta for reading!\n";
open (FASTA, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/TargetSeq.fa") || die "Can't open TargetSeq.fa for writing!\n";

my $TmpTit="";
my $TmpSeq="";
my $TmpTit_S="";
my $TmpTit_A="";
my $TmpSeq_Acomp;
my $TmpSeq_ARevcomp;

if ($Option eq "b") {              #B=both:sense and anti-sense strand

    while (<Inseq>) {

	  chomp $_; 
	 
      if (/^>(\S+)/){
	   
	    if ($TmpTit && $TmpSeq) {
	
			 $TmpTit_S=$TmpTit."_S";
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_S\n$TmpSeq\n";	
	         print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
		    }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	   }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	 
	}

       if ($TmpTit && $TmpSeq) {

	        my $TmpTit_S=$TmpTit."_S";
            my $TmpTit_A=$TmpTit."_A";

	        my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
     
            print FASTA ">$TmpTit_S\n$TmpSeq\n";
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n"; 
	
  }
}elsif ($Option eq "a") {           #Non-template: A: anti-sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	    if ($TmpTit && $TmpSeq) {
		  	  
		     $TmpTit_A=$TmpTit."_A";
		
			 $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
		     $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/;  	
		     $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
				
		     print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";	
		  }
	
	         $TmpTit=$1;	
	         $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	         $TmpSeq.=$1;
		
	  }
	}

       if ($TmpTit && $TmpSeq) {

            my $TmpTit_A=$TmpTit."_A";

            my $TmpSeq_Acomp=$TmpSeq;                        #Reverse-comp                                                
            $TmpSeq_Acomp =~ tr/atucgACGUT/TAAGCTGCAA/; 
            my $TmpSeq_ARevcomp = reverse($TmpSeq_Acomp); 
   
            print FASTA ">$TmpTit_A\n$TmpSeq_ARevcomp\n";
   
  }
}elsif ($Option eq "s" || $Option eq "") {     #Template: S: sense strand

	while (<Inseq>) {

		chomp $_;

	if (/^>(\S+)/){
	
	   if ($TmpTit && $TmpSeq) {
		  	  
		    $TmpTit_S=$TmpTit."_S";
				
		    print FASTA ">$TmpTit_S\n$TmpSeq\n";	
		  
	    }
	
	        $TmpTit=$1;	
	        $TmpSeq="";
	
	  }elsif(/(\w+)/){
	
	    	$TmpSeq.=$1;
		
	  } 
	}

      if ($TmpTit && $TmpSeq) {
	
		   my $TmpTit_S=$TmpTit."_S";

           print FASTA ">$TmpTit_S\n$TmpSeq\n";
    
  }
}
close(Inseq);
close(FASTA);

########################### Find CRISPR targets-single #################################
print "Step2: Find CRISPR targets.\n";
print  LOG "# Step2: Find CRISPR targets.\n";

open(IntPut, "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/TargetSeq.fa") || die "Can't open TargetSeq.fa for reading!\n";

open(OutPut, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_protospacer_single.txt") || die "Can't open report_protospacer_single.txt for writing!\n";
open(OutPut1, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa") || die "Can't open CRISPR.targets_single.fa for writing!\n";
#open(OutPut2, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/full_length_sgRNA.txt") || die "Can't open full_length_sgRNA.txt for writing!\n";
open(OutPut3, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_S.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";
open(OutPut4, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_A.txt") || die "Can't open CRISPR.targets_single.fa for writing!\n";

#print OutPut "Candidate sgRNA, Pattern: GGX18NGG, GX19NGG, X20NGG, GC% >=$GC %\n\n";
print OutPut "sgRID\t"."Start\t"."End\t"."CRISPR_target_sequence(5'-3')\t"."Length(nt)\t"."GC%\n";

my $ID=""; 
my $seq="";
my $dCas9handle = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCG";
my $terminator = "UUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU";
my $probe_id="";

while(<IntPut>) {
	chomp $_;

	if (/^>(\S+)/){
		analysis(); $ID=$1;

  }else{
		$seq=$_;
		
  }
}
analysis();
close OutPut;
close OutPut1;
close OutPut3;
close OutPut4;

############################# Find CRISPR targets-pairs #################################
open( PA, "<$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_A.txt" ) || die "can't open CRISPR.targets_A.txt!";
open( PS, "<$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_S.txt" ) || die "can't open CRISPR.targets_S.txt!";
open( PAIRS1, ">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_protospacer_pairs.xls" ) || die "can't open report_protospacer_pairs.xls!";
open( PAIRS2, ">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_pairs.fa" ) || die "can't open CRISPR.targets_pairs.fa!";

my $sgRID_A;
my $Start_A;
my $End_A;
my $target_seq_A ="";
my $Pattern_A;
my $GC_A;

my $sgRID_S;
my $Start_S;
my $End_S;
my $target_seq_S ="";
my $Pattern_S;
my $GC_S;

my $len_target_seq_A=0;
my $len_target_seq_S=0;

my $ID_A;
my $ID_S;

print PAIRS1 "\t\tPaired-gRNA\n";
print PAIRS1 "sgRID_S\ttarget_seq_S\tStart_S\tEnd_S\tGC%_S\t\tsgRID_A\ttarget_seq_A\tStart_A\tEnd_A\tGC%_A\tsgRNA_offset(bp)\n";

while ( <PA> ) {
	chomp $_; 
	(my $sgRID_A, my $Start_A, my $End_A, my $target_seq_A,	my $Pattern_A, my $GC_A, my $emp_A)=split/\t/, $_;

	$len_target_seq_A = length($target_seq_A);

	next if $len_target_seq_A ne ($truncat+3);

	$ID_A = $sgRID_A;
	$ID_A=~ s/_A_(\d+)//m;
	
	seek PS, 0, 0;

	while ( <PS> ) {
    chomp $_; 
		(my $sgRID_S, my $Start_S, my $End_S, my $target_seq_S,	my $Pattern_S, my $GC_S, my $emp_S)=split/\t/, $_;

		next if $target_seq_S eq "";
		$len_target_seq_S = length($target_seq_S);
		next if $len_target_seq_S ne ($truncat+3);               

		$ID_S = $sgRID_S;
		$ID_S=~ s/_S_(\d+)//m;
		
		my $offset_value = $Start_S -$End_A;

		if (($ID_A eq $ID_S) and ($offset_value > "$offset_s" and $offset_value < "$offset_e")) {      # -2 to 32 bp or 5 to 35 bp

			print PAIRS1 "$sgRID_A"."\t"."$target_seq_A"."\t"."$Start_A"."\t"."$End_A"."\t"."$GC_A"."\t<->\t";
			print PAIRS1 "$sgRID_S"."\t"."$target_seq_S"."\t"."$Start_S"."\t"."$End_S"."\t"."$GC_S"."\t"."$offset_value"."\n";
			print PAIRS2 ">"."$sgRID_A"."\n"."$target_seq_A"."\n";
			print PAIRS2 ">"."$sgRID_S"."\n"."$target_seq_S"."\n";

		}
  }
}
close(PA);
close(PS);
close(PAIRS1);
close(PAIRS2);

########################### Unique pairs sgR fasta seq ######################
my %seq;
my $title;
my $infile="$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_pairs.fa";

open (IN,"$infile") || die "can't open $infile!";
while (<IN>){
	$_=~s/\n//;
	$_=~s/\r//;
	if ($_=~/>/){
		$title=$_;
		$title=~s/>//;
	}
	else{
		$seq{$_}=$title;
	}
}
close IN;
#remove the abundant sequences
my @seq=keys (%seq);
my @uniqueseq;
my $find=0;
foreach (@seq){
	$find=0;
	my $seq=uc($_);
	foreach (@uniqueseq){
		if ($seq=~/$_/){
			$_=$seq;#replace with longer seq
			$find=1;
		}
		if ($_=~/$seq/){
			$find=1;
		}
	}
	if ($find==0){
		push @uniqueseq,$seq;
	}
}
#outout the final result
open (OUT,">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa");
foreach (@uniqueseq){
	print OUT ">$seq{$_}\n$_\n";
}
close OUT;

###################################### seqmap ###################################
print "Step3: Evaluate CRISPR potential off-target effect.\n";
print "Step3-1: Whole genome mapping.\n\n";
print  LOG "# Step3: Evaluate CRISPR potential off-target effect.\n";
print  LOG "# Step3-1: Whole genome mapping.\n";

if ($Seqmap_vesion eq "l" and $Type eq "s") {       #Linux-64

	my @args = ("./Seqmap/seqmap-1.0.12-linux-64","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "l" and $Type eq "p") {

	my @args = ("./Seqmap/seqmap-1.0.12-linux-64","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "w" and $Type eq "s") {   #windows 32, 64-bit

	my @args = ("./Seqmap/seqmap-1.0.12-windows.exe","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "w" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-windows.exe","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);

	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "m" and $Type eq "s") {    #macOSX-64

	my @args = ("./Seqmap/seqmap-1.0.12-mac-64","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
		die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "m" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-mac-64","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "a" and $Type eq "s") {    #macOSX-32

	my @args = ("./Seqmap/seqmap-1.0.12-mac","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "a" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-mac","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "u" and $Type eq "s") {    #Linux-32

	my @args = ("./Seqmap/seqmap-1.0.12-linux","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}elsif ($Seqmap_vesion eq "u" and $Type eq "p") {   

	my @args = ("./Seqmap/seqmap-1.0.12-linux","$Num_mismatch", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa", "$Genome", "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt", "/output_all_matches");
	system(@args);
	
	if($? == -1) {
	  die "system @args failed: $?";
	}

}

################################## count_sgR_Targetsite #########################
print "\nFinished mapping candidate CRISPR target sequences to genome: $Genome.\n";
print "\nStep3-2: Count the No.of potential off-target cleavage sites.\n";
print  LOG "# Finished mapping candidate CRISPR target sequences to genome: $Genome.\n";
print  LOG "# Step3-2: Count the No.of potential off-target cleavage sites.\n";

open(Input, "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt") || die "Can't open seqmap_output.txt for reading!\n";
open(Out, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_sgRid.txt") || die "Can't open report_sgRid.txt for writing!\n";

print Out "sgRID\t"."Total_No.OT\n";

my %ID;
my $len1=$truncat+2;
my $len2=$truncat+1;

while (<Input>) {

	chomp $_; 
	
	(my $trans_id, my $trans_coord, my $target_seq, my $probe_id, my $probe_seq, my $num_mismatch, my $strand)=split/\t/, $_;
	
	next if ($trans_id eq "trans_id" or (substr($target_seq,$len1,1) ne "G") or (substr($target_seq,$len2,1) eq "T") or (substr($target_seq,$len2,1) eq "C") ) ;

    while ($probe_id=~ /(\w[\w-]*)/g) {
		$ID{$1}++;
	}
}

foreach (sort {$ID{$a}<=>$ID{$b}}keys %ID) {
	
	next if $_=~ /probe_id/;
    
	print Out "$_\t$ID{$_}\n";

}
close Input;
close Out;

################################## format_seqmap #################################
print "Step3-3: Format and sort whole genome mapping result.\n";
print  LOG "# Step3-3: Format and sort whole genome mapping result.\n";

mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/A.Sort_OT_byID",0755)||die "can't create directory: $!" ;
#mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/D.Type_II_POT",0755)||die "can't create directory: $!" ;
mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/POT_seed.identity",0755)||die "can't create directory: $!" ;

my $trans_id;
my $trans_coord;
my $target_seq;
my $probe_id1;
my $probe_seq;
my $num_mismatch;
my $New_num_mismatch;
my $strand;

my $sgR_ID;

my $countsgR;

my $i=0;
my $j=0;

my $len3= $truncat+2; #19
my $len4= $truncat+1; #18
my $len5= $truncat-7; #10
my $len6= $truncat-12; #5

open(Seq, "<$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt" ) || die "can't open seqmap_output.txt!";
open(sgR, "<$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_sgRid.txt" ) || die "can't open report_sgRid.txt!";

open(Out3, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/search_OT.txt") || die "Can't open search_OT.txt for writing!\n";
print Out3 "\n\t\t------------------------------------ off-target analysis(OT:off-target) -----------------------------------\n\n";

while (<Seq>) {
	chomp $_; 

	(my $trans_id, my $trans_coord, my $target_seq, my $probe_id1, my $probe_seq, my $num_mismatch, my $strand)=split/\t/, $_;

	my $New_num_mismatch = $num_mismatch;
   
	seek sgR, 0, 0;
	
	while (<sgR>) {
		chomp $_; 
	
		(my $sgR_ID, my $countsgR)=split/\t/, $_;

		if ($probe_id1=~/$sgR_ID$/m) {
			$j++;
			    
			my @seq1 = split //, $probe_seq;
			my @seq2 = split //, $target_seq;
			
			my @Iden;
			$i++;
			
			for my $n (0..@seq1-1) {
				if   ( ($seq1 [$n] eq 'A' && $seq2 [$n] eq 'A') 
				    || ($seq1 [$n] eq 'T' && $seq2 [$n] eq 'T') 
				    || ($seq1 [$n] eq 'C' && $seq2 [$n] eq 'C') 
				    || ($seq1 [$n] eq 'G' && $seq2 [$n] eq 'G') ) {
				
				push @Iden, '-';
				}
				else {
				
					push @Iden, lc($seq2[$n]);
				
				}
			}  
   		#discard                                                        NG[N]                                  NTG                                  NCG
	    next if ($trans_id eq "trans_id" or (substr($target_seq,$len3,1) ne "G") or (substr($target_seq,$len4,1) eq "T") or (substr($target_seq,$len4,1) eq "C") ) ;
		
			open (Out2,">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/A.Sort_OT_byID/$sgR_ID.txt");
			
	   #                                                     NAG                      [N]-GG                                        
	    if($num_mismatch > 0 and (substr($target_seq,$len4,1) eq "A" or (substr($target_seq,$truncat,1) ne substr($probe_seq,$truncat,1)))){
	
	    	$New_num_mismatch=$num_mismatch-1;
	
				print Out2 "$probe_id1\t"."POT$i\t"."@Iden"."\t$New_num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
	    	#print Out3  "\t\t"."20           13 12      8 7           1 N G G\n";
	    	print Out3 "$probe_id1\t"."@seq1"."\n"."POT$i\t\t"."@seq2"."\n"."POT$i\t\t"."@Iden"."\t$New_num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
			}else {
	
				print Out2 "$probe_id1\t"."POT$i\t"."@Iden"."\t$num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
	    	#print Out3  "\t\t"."20           13 12      8 7           1 N G G\n";
	    	print Out3 "$probe_id1\t"."@seq1"."\n"."POT$i\t\t"."@seq2"."\n"."POT$i\t\t"."@Iden"."\t$num_mismatch"."M\t"."$trans_id\t"."$trans_coord\t"."$strand\n";
	
			}

			open (Out5,">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/POT_seed.identity/$sgR_ID.seed.identity.txt");
	    next if $strand =~/strand/;
	    #                                                12mer                                                       NAG                             [N]-GG 
	    if($num_mismatch > 0 and (substr($target_seq,$len6,12) eq substr($probe_seq,$len6,12)) and ((substr($target_seq,$len4,1) eq "A") or (substr($target_seq,$truncat,1) ne substr($probe_seq,$truncat,1)))){
	
	    	my $New3_num_mismatch=$num_mismatch-1;
		
	    	print Out5  "@Iden\t$New3_num_mismatch"."M\tPOT$i\t$trans_id\t$trans_coord\t$strand\t$probe_id1\n";
	
	    }elsif(substr($target_seq,$len6,12) eq substr($probe_seq,$len6,12)){
	
	    	print Out5  "@Iden\t$num_mismatch"."M\tPOT$i\t$trans_id\t$trans_coord\t$strand\t$probe_id1\n";
	 
			}
    }
  }
}
close(Seq);
close(sgR);
close(Out2);
close(Out3);

#################### count mismatch AND sort mapping by seed-ident.0-3Mismatch ################################
print "Step3-4: Sort mapping result by mismatch location.\n\n";
print  LOG "# Step3-4: Sort mapping result by mismatch location.\n\n";

mkdir("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/B.Sort_POT_byID",0755)||die "can't create directory: $!" ;

my $Dir = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/A.Sort_OT_byID" ;
my $Ext = "txt" ; #file type

opendir(DH, "$Dir") or die "Can't open: $!\n" ;

my @list = grep {/$Ext$/ && -f "$Dir/$_" } readdir(DH);

closedir(DH) ;
chdir($Dir) or die "Can't cd dir: $!\n" ;

my $file;

foreach my $file (@list){

	open(FH, "$file") || die "Can't open $file for reading!\n";
	
	$file =~ s/.txt//m;
	
	open(OUTi,">>$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/B.Sort_POT_byID/$file.POT.txt") || die "Can't open $file.POT.txt for writing!\n";
	
	my %ID;

	while (<FH>) {
 		chomp $_; 

  	(my  $sgRID, my $POT, my $SEQ, my $mismatch, my $Chr, my $Location, my $strand)=split/\t/, $_;

    while ($mismatch=~ /((\d)M|\w[\w-]*M)/g) {

	  	{
			$ID{$1}++;
			}
		}

		#next if ($mismatch eq "0M");
		#$truncat = 17, 18, 19, 20 nt sgRNA
		my $lentru17_7 = $truncat+3; #20
		my $lentru17_12= $truncat-7; #10
		
		my $lentru18_7 = $truncat+4; #22
		my $lentru18_12= $truncat-6; #12
		
		my $lentru19_7 = $truncat+5; #24
		my $lentru19_12= $truncat-5; #14
		
		my $lentru20_7 = $truncat+6; #26
		my $lentru20_12= $truncat-4; #16
		
		
		######                                                            7mer                                                 12mer
		#$truncat = 17
		if ($truncat ==17 and $mismatch =~/0M|1M|2M|3M/g and ((substr($SEQ,$lentru17_7,13) ne "- - - - - - -") or (substr($SEQ,$lentru17_12,23) ne "- - - - - - - - - - - -") )) {  
		
			print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."random_0_3M\t"."\n";
		
		}elsif ($truncat ==17 and (substr($SEQ,$lentru17_7,13) eq "- - - - - - -") and (substr($SEQ,$lentru17_12,23) ne "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."regionI_ident\t"."\n";
		
		}elsif ($truncat ==17 and(substr($SEQ,$lentru17_12,23) eq "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."seed_ident\t"."\n";
		
		 #$truncat = 18
		}elsif ($truncat ==18 and $mismatch =~/0M|1M|2M|3M/g and ((substr($SEQ,$lentru18_7,13) ne "- - - - - - -") or (substr($SEQ,$lentru18_12,23) ne "- - - - - - - - - - - -") )) {  
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."random_0_3M\t"."\n";
		
		}elsif ($truncat ==18 and (substr($SEQ,$lentru18_7,13) eq "- - - - - - -") and (substr($SEQ,$lentru18_12,23) ne "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."regionI_ident\t"."\n";
		
		}elsif ($truncat ==18 and(substr($SEQ,$lentru18_12,23) eq "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."seed_ident\t"."\n";
		
		 #$truncat = 19
		}elsif ($truncat ==19 and $mismatch =~/0M|1M|2M|3M/g and ((substr($SEQ,$lentru19_7,13) ne "- - - - - - -") or (substr($SEQ,$lentru19_12,23) ne "- - - - - - - - - - - -") )) {  
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."random_0_3M\t"."\n";
		
		}elsif ($truncat ==19 and (substr($SEQ,$lentru19_7,13) eq "- - - - - - -") and (substr($SEQ,$lentru19_12,23) ne "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."regionI_ident\t"."\n";
		
		}elsif ($truncat ==19 and(substr($SEQ,$lentru19_12,23) eq "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."seed_ident\t"."\n";
		
		 #$truncat = 20
		}elsif ($truncat ==20 and $mismatch =~/0M|1M|2M|3M/g and ((substr($SEQ,$lentru20_7,13) ne "- - - - - - -") or (substr($SEQ,$lentru20_12,23) ne "- - - - - - - - - - - -") )) {  
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."random_0_3M\t"."\n";
		
		}elsif ($truncat ==20 and (substr($SEQ,$lentru20_7,13) eq "- - - - - - -") and (substr($SEQ,$lentru20_12,23) ne "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."regionI_ident\t"."\n";
		
		}elsif ($truncat ==20 and(substr($SEQ,$lentru20_12,23) eq "- - - - - - - - - - - -")) {
		
		  print OUTi "$SEQ\t"."$mismatch\t"."$sgRID"."_$POT\t"."$Chr\t"."$Location\t"."$strand\t"."seed_ident\t"."\n";
		}
  }

######
    foreach (sort {$ID{$a} <=> $ID{$b}}keys %ID) {
 
			#print OUTM "$_"."\t"."$ID{$_}"."\t";

    }
    #print OUTM "\n";

}
#close(FH) ;
#close (OUTM);


################################### subroutine ###################################
sub analysis {  
	if($ID eq "" || $seq eq "") {
	  return();
	}
	#my $seq =~ s/\n//g;
	my $SEQRNA  = uc($seq);
	my $SEQDNA  = $SEQRNA;
	$SEQDNA  =~ tr/U/T/;

	my $len = length($seq);

	my $idx_s;
	my $i = 0;

	for($idx_s=-1; $idx_s<length($seq); $idx_s++) {

		my $lentruncat = $truncat+3;
		my $lentruncat2 = $truncat+1;

		if (substr($SEQDNA, $idx_s+1, $lentruncat) =~ /[ATCG]{$lentruncat2}GG/) {    #X?-NGG  
			my $sgRidx_s = $idx_s+2;                                  #For sense strand
			my $sgR = $&;
			#print $sgR."\n";
			my $SgR1 = substr($sgR,0,-3); 
			#print $SgR1."\n";
			my $sgRAT = &AT($sgR);
			#my $sgRGC = &GC($sgR);
			my $sgRGC = &GC($SgR1);

			#my $SgRsd = substr($sgR,8,12);                           #5-X12-NGG-3
			#my $SgRseed ="$SgRsd"."NGG";
			my $SgRNt = substr($sgR,0,$truncat);                      #truncat
			my $sgRlen = length($sgR);
			my $sgRidx_e = $sgRidx_s+$sgRlen-1;
			$i++;
			
			my $sgRidx_A_s =$len-$sgRidx_e+1;                         #For anti-sense strand
			my $sgRidx_A_e =$len-$sgRidx_s+1;
			
			my $sgrna = $SgRNt;
			$sgrna =~ tr/T/U/;
			my $sgRNA = $sgrna.$dCas9handle.$terminator;

				if (!($sgR =~ /T{4,18}/g)){
					#if (!($sgR =~ /A{5,21}|T{4,21}|C{6,21}|G{6,21}|(AT){6,10}|(AC){6,10}|(AG){6,10}|(TA){6,10}|(TC){6,10}|(TG){6,10}|(CA){6,10}|(CT){6,10}|(CG){6,10}|(GA){6,10}|(GT){6,10}|GC{6,10}|(AAT){5,7}|(AAC){5,7}|(AAG){5,7}|(ATA){5,7}|(ATT){5,7}|(ATC){5,7}|(ATG){5,7}|(ACA){5,7}|(ACT){5,7}|(ACC){5,7}|(ACG){5,7}|(AGA){5,7}|(AGT){5,7}|(AGC){5,7}|(AGG){5,7}|(TAA){5,7}|(TAT){5,7}|(TAC){5,7}|(TAG){5,7}|(TTA){5,7}|(TTC){5,7}|(TTG){5,7}|(TCA){5,7}|(TCT){5,7}|(TCC){5,7}|(TCG){5,7}|(TGA){5,7}|(TGT){5,7}|(TGC){5,7}|(TGG){5,7}|(CAA){5,7}|(CAT){5,7}|(CAC){5,7}|(CAG){5,7}|(CTA){5,7}|(CTT){5,7}|(CTC){5,7}|(CTG){5,7}|(CCA){5,7}|(CCT){5,7}|(CCG){5,7}|(CGA){5,7}|(CGT){5,7}|(CGC){5,7}|(CGG){5,7}|(GAA){5,7}|(GAT){5,7}|(GAC){5,7}|(GAG){5,7}|(GTA){5,7}|(GTT){5,7}|(GTC){5,7}|(GTG){5,7}|(GCA){5,7}|(GCT){5,7}|(GCC){5,7}|(GCG){5,7}|(GGA){5,7}|(GGT){5,7}|(GGC){5,7}/g) ){  

					if ($ID=~ /._S$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m)) {    #For sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						print OutPut3 "$ID\_$i\t"."$sgRidx_s\t"."$sgRidx_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";

					}elsif ($ID=~ /._A$/ and ($sgRGC >= $GC_l and $sgRGC < $GC_m) ) {   #For anti-sense strand
           
						print OutPut "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						print OutPut4 "$ID\_$i\t"."$sgRidx_A_s\t"."$sgRidx_A_e\t"."$sgR\t"."$sgRlen\t"."$sgRGC %\n";
						
						print OutPut1 ">$ID\_$i\n"."$sgR\n";
						#print OutPut2 ">$ID\_$i\n"."$sgRNA"."\n";

    			}
				}
		}
	}
}

########count_ot#########
my $Dir3_ot = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/A.Sort_OT_byID" ;
my $Ext3_ot = "txt" ; #file type

my ($OTname_ot, $OTID_ot, $pattern_ot, $mismatch_ot,$chrosome_ot,$location_ot,$strand_ot);

opendir(DH, "$Dir3_ot") or die "Can't open: $!\n" ;
my @list_ot = grep {/$Ext3_ot$/ && -f "$Dir3_ot/$_" } readdir(DH) ;
closedir(DH) ;
chdir($Dir3_ot) or die "Can't cd dir: $!\n" ;

my $file_ot;

open(OutPut, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_ot.txt") || die "Can't open report_count_total_ot.txt for writing!\n";

print OutPut "ID\t0M(on-/off-)\t1M\t2M\t3M\t4M\t5M\tTotal OT number\n";

foreach my $file_ot (@list_ot){
    open(FH, "$file_ot") or die "Can't open: $!\n" ;
	$file_ot =~ s/.txt//m;

my $j=0;
my $e=0;
my $count_0M=0;
my $count_1M=0;
my $count_2M=0;
my $count_3M=0;
my $count_4M=0;
my $count_5M=0;
my $line1_ot;
my $count_error;

while($line1_ot=<FH>){
    $j++;
	($OTname_ot, $OTID_ot, $pattern_ot, $mismatch_ot,$chrosome_ot,$location_ot,$strand_ot)=split(/\t/, $line1_ot);
	#pig_NR1H4_cds_A_4	POT264	- - - - - - - t - - - - - - - t g - - - g - -	3M	10	9370885	+
	if ($mismatch_ot eq '0M'){
	 ++$count_0M;
	} elsif ($mismatch_ot eq '1M'){
	 ++$count_1M;
	} elsif ($mismatch_ot eq '2M'){
	 ++$count_2M;
	} elsif ($mismatch_ot eq '3M'){
	 ++$count_3M;
	} elsif ($mismatch_ot eq '4M'){
	 ++$count_4M;
	} elsif ($mismatch_ot eq '5M'){
	 ++$count_5M;
	} else {
	 print OutPut "Error-Strange mismatch numbers:$mismatch_ot\n";
	 ++$count_error;
	}	
}

$e=$j-1;

    print OutPut $file_ot."\t";
    print OutPut "$count_0M\t";
	print OutPut "$count_1M\t";
	print OutPut "$count_2M\t";
	print OutPut "$count_3M\t";
    print OutPut "$count_4M\t";
	print OutPut "$count_5M\t";

if ($count_0M == 0) {

	print OutPut "$j\n";

}else {

    print OutPut "$e\n";
}
}

close OutPut;


########count_pot#########
my $Dir3 = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/POT_seed.identity" ;
my $Ext3 = "txt" ; #file type

my $pattern;
my $mismatch;
my $name;
my $chrosome;
my $location;
my $starnd;
my $cds_name;

opendir(DH, "$Dir3") or die "Can't open: $!\n" ;
my @list_pot = grep {/$Ext3$/ && -f "$Dir3/$_" } readdir(DH) ;
closedir(DH) ;
chdir($Dir3) or die "Can't cd dir: $!\n" ;

my $file_pot;
open(OutPut, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_pot.txt") || die "Can't open report_count_total_pot.txt for writing!\n";

print OutPut "ID_seed.identity\t0M(on-/off-)\t1M\t2M\t3M\t4M\t5M\tTotal POT number\n";

foreach my $file_pot (@list_pot){
    open(FH, "$file_pot") or die "Can't open: $!\n" ;
	$file_pot =~ s/.seed.identity.txt//m;

my $j=0;
my $e=0;
my $count_0M=0;
my $count_1M=0;
my $count_2M=0;
my $count_3M=0;
my $count_4M=0;
my $count_5M=0;
my $line1;
my $count_error;

while($line1=<FH>){
    $j++;
	my($pattern, $mismatch, $name, $chrosome, $location, $starnd, $cds_name)=split(/\t/, $line1);
	if ($mismatch eq '0M'){
	 ++$count_0M;
	} elsif ($mismatch eq '1M'){
	 ++$count_1M;
	} elsif ($mismatch eq '2M'){
	 ++$count_2M;
	} elsif ($mismatch eq '3M'){
	 ++$count_3M;
	} elsif ($mismatch eq '4M'){
	 ++$count_4M;
	} elsif ($mismatch eq '5M'){
	 ++$count_5M;
	} else {
	 print OutPut "Error-Strange mismatch numbers:$mismatch\n";
	 ++$count_error;
	}	
}

$e=$j-1;
    print OutPut $file_pot."\t";
    print OutPut "$count_0M\t";
	print OutPut "$count_1M\t";
	print OutPut "$count_2M\t";
	print OutPut "$count_3M\t";
    print OutPut "$count_4M\t";
	print OutPut "$count_5M\t";

if ($count_0M == 0) {

    print OutPut "$j\n";

}else {

    print OutPut "$e\n";
}
}
close OutPut;


#########
#########
my $File1_ot = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_protospacer_single.txt";
my $File2_ot = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_ot.txt";

my ($sgRID_ot,	$Start_ot, $End_ot, $CRISPR_target_sequence_ot,$Length_ot,$GC_ot);
my ($ID_ot, $M0_ot, $M1_ot,$M2_ot, $M3_ot, $M4_ot,$M5_ot,$OT_ot);

open( PA, "$File1_ot" ) || die "can't open $File1_ot!";
open( PS, "$File2_ot" ) || die "can't open $File2_ot!";
open( PAIRS1, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/single_and_ot.txt" ) || die "can't open single_and_ot.txt!";

while ( <PA> ) {
    chomp();
($sgRID_ot,$Start_ot, $End_ot, $CRISPR_target_sequence_ot,$Length_ot,$GC_ot)=split/\t/, $_;

    seek PS, 0, 0;

    while (  <PS> ) {
    chomp();
   ($ID_ot, $M0_ot, $M1_ot,$M2_ot, $M3_ot, $M4_ot,$M5_ot,$OT_ot)=split/\t/, $_;

if ($ID_ot eq $sgRID_ot) {

print PAIRS1 $sgRID_ot."\t".$Start_ot."\t".$End_ot."\t".$CRISPR_target_sequence_ot."\t".$Length_ot."\t".$GC_ot."\t".$ID_ot."\t".$M0_ot."\t".$M1_ot."\t".$M2_ot."\t".$M3_ot."\t".$M4_ot."\t".$M5_ot."\t".$OT_ot."\n";

}
}
}
close(PA);
close(PS);
close(PAIRS1);

#########
my $File1 = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/single_and_ot.txt";
my $File2 = "$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_pot.txt";

my ($sgRID_POT,	$Start_POT, $End_POT, $CRISPR_target_sequence_POT,$Length_POT,$GC_POT,$ID_POT, $M0_POT, $M1_POT,$M2_POT, $M3_POT, $M4_POT,$M5_POT,$OT_POT);
my ($ID_POT2, $M0_POT2, $M1_POT2,$M2_POT2, $M3_POT2, $M4_POT2,$M5_POT2,$OT_POT2,$Risk_evaluation);


open( PA, "$File1" ) || die "can't open $File1!";
open( PS, "$File2" ) || die "can't open $File2!";
open( PAIRS1, ">$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/sgRNAcas9_report.xls" ) || die "can't open sgRNAcas9_report.xls!";

print PAIRS1 "sgRID\t"."Start\t"."End\t"."Protospacer_sequence+PAM(5'-3')\t"."Length(nt)\t"."GC%_of_Protospacer\t"."Protospacer+PAM(OT)\t"."0M(on-/off-)\t"."1M\t"."2M\t"."3M\t"."4M\t"."5M\t"."Total_No.of_OT\t"."Seed_12+PAM(POT)\t"."0M(on-/off-)\t"."1M\t"."2M\t"."3M\t"."4M\t"."5M\t"."Total_No.of_POT\t"."Risk_evaluation\n";

while ( <PA> ) {
    chomp();
($sgRID_POT, $Start_POT, $End_POT, $CRISPR_target_sequence_POT,$Length_POT,$GC_POT,$ID_POT, $M0_POT, $M1_POT,$M2_POT, $M3_POT, $M4_POT,$M5_POT,$OT_POT)=split/\t/, $_;

#print "$sgRID\n";

    seek PS, 0, 0;

    while ( <PS> ) {
    chomp();
   ($ID_POT2, $M0_POT2, $M1_POT2,$M2_POT2, $M3_POT2, $M4_POT2,$M5_POT2,$OT_POT2)=split/\t/, $_;

#print "$ID\n";

if ($ID_POT2 eq $sgRID_POT) {

if ($M0_POT == 0 ) {

$Risk_evaluation = "Discard";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t".""."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}elsif ($M0_POT >=2 ) {

$Risk_evaluation = "Repeat_sites_or_bad?";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}elsif ($M1_POT ==0 and $M1_POT2 ==0 and $M2_POT ==0 and $M2_POT2 ==0 and $M3_POT2 ==0 and $M4_POT2 ==0 and $M5_POT2 ==0 ){

$Risk_evaluation = "Best";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}elsif ($M1_POT ==0 and $M1_POT2 ==0 and $M2_POT ==0 and $M2_POT2 ==0 and $M3_POT2 >=0 and $M4_POT2 >=0 and $M5_POT2 >=0 ){

$Risk_evaluation = "Low_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}elsif ($M1_POT>=1 or $M1_POT2 >=1 ) {

$Risk_evaluation = "High_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}elsif ($M2_POT>=1 or $M2_POT2 >=1 ){

$Risk_evaluation = "Moderate_risk";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";

}else{

$Risk_evaluation = "Unclassified";

print PAIRS1 $sgRID_POT."\t".$Start_POT."\t".$End_POT."\t".$CRISPR_target_sequence_POT."\t".$Length_POT."\t".$GC_POT."\t"."#"."\t".$M0_POT."\t".$M1_POT."\t".$M2_POT."\t".$M3_POT."\t".$M4_POT."\t".$M5_POT."\t".$OT_POT."\t"."#"."\t".$M0_POT2."\t".$M1_POT2."\t".$M2_POT2."\t".$M3_POT2."\t".$M4_POT2."\t".$M5_POT2."\t".$OT_POT2."\t".$Risk_evaluation."\n";
}
}
}
}
close(PA);
close(PS);
close(PAIRS1);

#######subroutine to calculate GC% content
sub GC { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /G|C/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}


#######subroutine to calculate AT% content
sub AT { 

    my $seq2 = shift @_; 

    $seq2 = uc($seq2);

    my $seq2DNA;
    $seq2DNA = $seq2;
    $seq2DNA =~ tr/U/T/;

    my $seq2length = length($seq2DNA);

    my $count = 0;
    for (my $i = 0; $i < $seq2length; $i++) {
			my $sub = substr($seq2DNA,$i,1);
			if ($sub =~ /A|T/i) {
	    	$count++;
			}
    }
    my $gc = sprintf("%.1f",$count * 100 /$seq2length);
    return $gc;
}

close IntPut;
close OutPut;
#exit;


#delete Temporary files.
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/seqmap_output.txt")||die "Can't delete seqmap_output.txt file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/search_OT.txt")||die "Can't delete search_OT.txt file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/unique_pairs.fa")||die "Can't delete unique_pairs.fa file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_A.txt")||die "Can't delete CRISPR.targets_A.txt file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_S.txt")||die "Can't delete CRISPR.targets_S.txt file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_pairs.fa")||die "Can't delete CRISPR.targets_pairs.fa file";
#unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/CRISPR.targets_single.fa")||die "Can't delete CRISPR.targets_single.fa file";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_sgRid.txt")||die "Can't delete report_sgRid.txt file";
#unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_protospacer_single.txt")||die "Can't delete report_protospacer_single.txt";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_ot.txt")||die "Can't delete report_count_total_ot.txt";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/report_count_total_pot.txt")||die "Can't delete report_count_total_pot.txt";
unlink ("$dir/sgRNAcas9.report_$truncat.$Option.$Inputfile_Fasta/single_and_ot.txt")||die "Can't delete single_and_ot.txt";


######################################## Running Time #################################################
my $oo = time() -  $oldtime;
printf  "\nTotal time consumption is $oo second.\n";

printf LOG "# Total time consumption is $oo seconds." ."\n"; 
print "Your job is done, open sgRNAcas9.report_$truncat.$Option directory to check result.". "\n";
print  LOG "# Your job is done, open sgRNAcas9.report_$truncat.$Option directory to check result.". "\n";

print  LOG "################################# END ###########################################"."\n\n";



