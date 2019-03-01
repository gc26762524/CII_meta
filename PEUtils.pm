package PEUtils;
use strict;
use warnings;
use File::Basename;

############################################################################################################################################
#
#	PEUtils.pm : This module contains subtroutines for paired end files and scripts
#
#############################################################################################################################################

sub get_fastq_file_mate_name {

    my $fastq = $_[0];
    my $fastq_mate = $fastq;
    $fastq_mate =~ s/\.R1/\.R2/g;
    return $fastq_mate;
}

#Return absolute path to Stitched directory
sub get_stitched_dir {

    my $fastq = $_[0];
    my ($fileName, $dir) = fileparse($fastq);
    my $fileDirName = basename($dir);
    $dir =~ s/$fileDirName\///;
    my $stitched_dir = "$dir"."Stitched"; #dir already includes a slash in the name
    #print "Stitched directory location: $stitched_dir\n";
    return $stitched_dir;

}

#Return Sample base name for current sample
sub get_sample_base_name_4stitching {

    my $fastq = $_[0];
    my $stitching_software = $_[1];

    my ($fileName, $fileDir) = fileparse($fastq);

    my $baseName = $fileName;
    $baseName =~ s/\.R1\.f[qa](stq)*/_${stitching_software}/g;

    my @suffixlist = (".fq",".fq.gz",".fastq",".fastq.gz",".gz");
    my $baseName = basename($baseName, @suffixlist);  
    
    return $baseName;

}

#Return Absolute File Names used for current sample in stitching step
sub get_out_fileNames_4stitching {

    my $fastq = $_[0];
    my $stitching_software = $_[1];

    my $stitched_dir = get_stitched_dir($fastq);
    print ("Stitched Dir: $stitched_dir\n");

    my $fileName = fileparse($fastq);

    my $outJoin = "$stitched_dir/".$fileName;
    $outJoin =~ s/\.R1\.f[qa](stq)*/_${stitching_software}.join.fq/g;

    my $outS1=$outJoin;
    $outS1=~s/join/S1/g;

    my $outS2=$outS1;
    $outS2=~s/S1/S2/g;

    my $outLogLength=$outJoin;
    $outLogLength=~s/fq/log/g;

    my $outLogResults=$outJoin;
    $outLogResults=~s/join\.fq/log/g;

    return ($outJoin, $outS1, $outS2, $outLogLength, $outLogResults);

}

1;
