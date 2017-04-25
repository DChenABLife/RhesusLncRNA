#!/usr/bin/perl -w
# 

my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your programme，you must write the detailed time、discriptions、parameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions(\%opts,"mrna=s","o=s", "d=s","l=s","mirna=s","c=s","h");

if(!defined($opts{mrna}) || !defined($opts{o}) || !defined($opts{d}))
{
	print <<"	Usage End.";

	Description:This programme is used for ~
		
		Version: $ver

	Usage:perl $0 -i

		-mrna      mrnafile                          must be given
		-o      outfile                         must be given
		-d      outdirectory                    must be given 
		-l      cut off file per lines[10]    choice
		-mirna      mirnafile              must be given
		-c      the correlation method,pearson:0,spearman:1 default is 0

	Usage End.

	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################
# >Locus_1_Transcript_1/1_Confidence_1.000_Length_559#
my $mrna=$opts{mrna};
my $out=$opts{o};
my $outdir=$opts{d};
my $line=$opts{l} || 100;
my $mirna=$opts{mirna};
my $method=$opts{c} || "0";

# $out=&AbsolutePath("file",$out);
$outdir = defined $outdir ? $outdir : "./";
`mkdir -p $outdir ` unless(-d $outdir);
$outdir=&AbsolutePath("dir",$outdir);

my $subfile = "";
my $linenum = 0;
my $subcount = 0;
open (IN,"$mrna") || die "Can't open $mrna\n";	
my $title = <IN>;
chomp($title);
close IN;
# open (TITLE,"title");
# print TITLE $title;
# close TITLE;


my $outcome =`wc -l $mrna `;
print "$outcome\n";
my $linecount = (split/\s+/,$outcome)[0];
my $a = length(int($linecount/$line));
my $out_name = basename($mrna);
`split -l $line -d -a $a $mrna $outdir/$out_name\_subset `;			#split the fa file into subset fa file
my $raw_postfix = basename($mrna)."_subset";
chdir($outdir);
my @raw;
&load_raw_data($outdir,\@raw,$raw_postfix);
print join("\n",@raw),"\n";
die("not found the raw data\n") if (!@raw);

open SH,">cor.sh";

for (my $i=0;$i<@raw ;$i++) {
	`sed -i "1i\\$title" $outdir/$raw[$i]` if $i!=0;
	print SH "cd $outdir && R --slave < $bin/cor.r --args $outdir/$raw[$i] $mirna cor_result.txt.$i p_result.txt.$i $method ./ && \n";
}
close SH;

`perl /public/bin/qsub-sge.pl --queue all.q --resource vf=10.0G --maxproc 600 cor.sh`;

`cat cor_result.txt.* > $out\_cor_result.tmp && rm -rf cor_result.txt.*`;
`cat p_result.txt.* > $out\_p_result.tmp && rm -rf p_result.txt.*`;

open SH,">cor_h.sh";

print SH "cat $out\_cor_result.tmp | perl -ne 'BEGIN{\$flag=0;}chomp;if(\/\^mRNA\/){if(\$flag==0){\$flag++;}else{next;}}print \$_,\"\\n\";' > $out\_cor_result.txt \n";
print SH "cat $out\_p_result.tmp | perl -ne 'BEGIN{\$flag=0;}chomp;if(\/\^mRNA\/){if(\$flag==0){\$flag++;}else{next;}}print \$_,\"\\n\";' > $out\_p_result.txt \n";
close SH;

`sh cor_h.sh`;

`rm -rf $out_name\_subset*`;

die ;

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub load_raw_data{
	my ($INDIR,$filenames_ref,$RAW_POSTFIX)=@_;
	#print "($INDIR,$filenames_ref,$RAW_POSTFIX)\n";die;
	opendir(DIR,$INDIR) or die "Can't open $INDIR: $!";
	my $tmp;
	while ($tmp=readdir(DIR)) {
		chomp $tmp;
		next if ($tmp!~/$RAW_POSTFIX/) ;
		push @{$filenames_ref},$tmp;
	}
	@{$filenames_ref} = sort @{$filenames_ref};
	close(DIR);
}
sub AbsolutePath{		#获取指定目录或文件的决定路径
	my ($type,$input) = @_;
	my $return;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	}
	elsif($type eq 'file'){
		my $pwd = `pwd`;
		chomp $pwd;
		my $dir=dirname($input);
		my $file=basename($input);
		chdir($dir);
		$return = `pwd`;
		chomp $return;
		$return .="\/".$file;
		chdir($pwd);
	}
	return $return;
}
