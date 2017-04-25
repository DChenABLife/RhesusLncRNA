#!/usr/bin/perl -w
# 
# Copyright (c)   ABLife
# Programmer:     CHENGCHAO <chaocheng@ablife.cc>
# 

=head1 Name

filter_target.pl

=head1 Description

=head1 Version

=head1 Usage
  
  -t           pure targets file                       input data directery, must be given;
  -c           correlation file                    output data directery,default is ./
  -pvalue      pvalue                             /users/xux/arabidopsis/TAIR10/Genome_sequence/TAIR10
  -cf          correlation filter                               /users/xux/test/Mus_musculus/bowtie_index/mm9
  -pf          pvalue filter                       
  -o           outfile                             


=head1 Exmple

  perl $0 -i /data2/arabidopsis/RNA-seq/repeat2/ -o _0921 -db /users/xux/arabidopsis/TAIR10/Genome_sequence/TAIR10 -p _0921 -gff /data2/arabidopsis/TAIR10/GFF/TAIR10_GFF3_genes.gff -r '.fq' -kgXref /data2/arabidopsis/TAIR10/function_annotation/gene_aliases.20101027 &


=cut
my $ver = "1.0.0";


use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use FileHandle;

my $buf = "-";
my %opts;
GetOptions( \%opts,"c=s","pvalue=s","cf=s","pf=s","o=s", "h" );

if (   !defined( $opts{c} )
	|| !defined( $opts{pvalue} )
	|| !defined( $opts{cf} )
	|| !defined( $opts{pf} )
	|| !defined( $opts{o} )
	|| defined( $opts{h} ) )
{
	print <<"	Usage End.";

		Version:$ver

	Usage:perl $0

        -c           correlation file                   
        -pvalue      pvalue
        -cf          correlation filter
        -pf          pvalue filter
        -o           outfile

	Usage End.

	exit;
}

#############Time_start#############
my $start_time = time();

my $Time_Start;
$Time_Start = sub_format_datetime( localtime( time() ) );
print "\nStart Time :[$Time_Start]\n\n";
####################################

#my $target_file  = $opts{t};
my $cor_file  = $opts{c};
my $pvalue_file  = $opts{pvalue};
my $cor_filter  = $opts{cf};
my $pvalue_filter  = $opts{pf};
my $out_file = $opts{o};

##读取target文件
#miRNA gene miranda_score targetscan_score
my %targets_hash = ();
#open TARGETS,"$target_file" || die;
#while(<TARGETS>){
#	chomp;
#	my @line = split(/\t/);
#	$targets_hash{$line[0]}->{$line[1]}="$line[2]\t$line[3]";
#}
#close TARGETS;



##读取correlation文件
#mRNA/miRNA(pvalue)      hsa-let-7a      hsa-let-7b      hsa-let-7c      hsa-let-7d
#AFFX-BioB-3     0.250165056525374       0.187708411987572       0.402625489886401       0.513094
my %cor_hash = ();
open COR,"$cor_file" || die;
my $head = <COR>;
chomp($head);
my @mirna=split(/\t/,$head);
while(<COR>){
	chomp;
	my @line=split(/\t/);
	my $gene=$line[0];
	for(my $i=1;$i<@line;$i++){
		next if($line[$i] eq "NA");				#skip the NA element
		next if(abs($line[$i])<=$cor_filter);
		next if($mirna[$i] eq $gene);		#skip the same gene.
		$cor_hash{$mirna[$i]}->{$gene}=$line[$i];
	}
}
close COR;
print "READ COR Done!\n";

##读取pvalue文件
#
my %pvalue_hash=();
open PVALUE,"$pvalue_file" || die;
$head = <PVALUE>;
chomp($head);
@mirna=split(/\t/,$head);
while(<PVALUE>){
	chomp;
	my @line=split(/\t/);
	my $gene=$line[0];
	for(my $i=1;$i<@line;$i++){
		next if($line[$i] eq "NA");				#skip the NA element
		next if($line[$i]>$pvalue_filter);
		next if($mirna[$i] eq $gene);		#skip the same gene.
		$pvalue_hash{$mirna[$i]}->{$gene}=$line[$i];
	}
}
close PVALUE;
print "Read PVALUE Done!\n";

open OUT,">$out_file";
print OUT "#lncRNA\tgene\tcorrelation\tpvalue\n";

##筛选出符合条件的target并输出
#miRNA gene cor pvalue miranda_score targetscan_score
foreach my $miRNA (keys %cor_hash){
	foreach my $gene(keys %{$cor_hash{$miRNA}}){
#next if (not defined $cor_hash{$miRNA}{$gene});
		next if(not defined $pvalue_hash{$miRNA}{$gene});
#		if(abs($cor_hash{$miRNA}{$gene})>=$cor_filter && $pvalue_hash{$miRNA}{$gene}<$pvalue_filter){
			print OUT $miRNA,"\t",$gene,"\t",$cor_hash{$miRNA}{$gene},"\t",$pvalue_hash{$miRNA}{$gene},"\n";
#		}
	}
}

############Time_end#############
my $Time_End;
$Time_End = sub_format_datetime( localtime( time() ) );
print "\nEnd Time :[$Time_End]\n\n";

my $time_used = time() - $start_time;
my $h = $time_used/3600;
my $m = $time_used%3600/60;
my $s = $time_used%3600%60;
printf("\nAll Time used : %d hours\, %d minutes\, %d seconds\n\n",$h,$m,$s);


#######Sub_format_datetime#######
sub sub_format_datetime {    #Time calculation subroutine
	my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = @_;
	$wday = $yday = $isdst = 0;
	sprintf(
		"%4d-%02d-%02d %02d:%02d:%02d",
		$year + 1900,
		$mon + 1, $day, $hour, $min, $sec
	);
}
