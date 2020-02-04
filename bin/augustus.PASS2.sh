#!/bin/bash

prefix=$1
fasta=$2
proteins=$3
aed=$4
local_config=$5
training_params=$6

AUGUSTUS_CONFIG_PATH=$local_config

basedir=$(pwd)
#prepare training parameters from BUSCO
if [ ! -z $training_params ]
then
	echo -e "[$(date)]\tPreparing training set from previous Augustus run"
	#get local copy of Augustus parameters from previous training round
	cp -fr $training_params $local_config/species/$prefix

	#rename files to current prefix
	cd $local_config/species/$prefix
	base=$(ls *weightmatrix.txt | sed 's/_weightmatrix.txt//')
	#rename files
	for file in $(ls -1); do new=$(echo -e "$file" | sed "s/$base/$prefix/g"); mv $file $new; done
	#rename the files cited within certain HMM configuration files
	sed -i "s/$base/$prefix/g" $prefix\_parameters.cfg
	sed -i "s/$base/$prefix/g" $prefix\_parameters.cfg.orig1

	cd $basedir
fi

if [ ! -z $aed ]
then
	echo -e "[$(date)]\tFiltering proteins with AED > $aed"
	#extract only proteins with AED < x
	cat <(echo -e "$aed") <(cat $proteins | perl -ne 'chomp; if ($_ =~ /^>/){print "\n$_\n"}else{print "$_"}' | grep -v "^$") | \
perl -ne 'chomp; if ($. == 1){$AED = $_}else{$h=$_; $s=<>; @a=split(" "); $a[2] =~ s/AED://; if ($a[2] < $AED){print "$h\n$s"}}' | sed 's/ .*//' > $prefix.AED-st$aed.maker.proteins.fasta
	proteins=$prefix.AED-st$aed.maker.proteins.fasta
fi

echo -e "[$(date)]\tRun AUGUSTUS"
autoAug.pl --genome=$fasta --species=$prefix --trainingset=$proteins --singleCPU -v --useexisting
retVal=$?

if [ ! $retVal -eq 0 ]
then
	>&2 echo "Augustus ended in an error"
	exit $retVal
fi

#copy the training set that was produced
cp -rf $local_config/species/$prefix .

echo -e "[$(date)]Reformatting to $(pwd)/autoAug/autoAugPred_abinitio/predictions/augustus.gff to GFF3 -> $(pwd)/augustus.gff3"
cat autoAug/results/predictions/augustus.abinitio.gff | perl -ne 'chomp; @a=split(/\t/); if ($a[2] eq 'gene'){$id=$a[-1]; $a[-1] =~ s/^/ID=/; print join("\t", @a)."\n"}else{if ($_ =~ /;$/){print "$_ Parent=$id\n"}else{print "$_; Parent=$id\n"}}' | sed 's/; /;/g' | sed 's/ /=/g' > augustus.gff3

#cat $(pwd)/autoAug/autoAugPred_abinitio/predictions/augustus.gff | perl -ne 'chomp; @a=split(/\t/); if ($a[2] eq 'gene'){$id=$a[-1]; $a[-1] =~ s/^/ID=/; print join("\t", @a)."\n"}else{if ($_ =~ /;$/){print "$_ Parent=$id\n"}else{print "$_; Parent=$id\n"}}' | sed 's/; /;/g' | sed 's/ /=/g' > $(pwd)/autoAug/autoAugPred_abinitio/predictions/augustus.gff3 

