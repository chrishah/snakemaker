configfile: "data/config.yaml"

import pandas as pd
import os
import glob
from math import ceil
from pathlib import Path

localrules: initiate, all

n=int(config["split_batch_length"])

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]

dic = {'sample': [], 'unit': []}

#function that partitions the fasta file
def parition_by_length(fasta, max_length=1000000, pr=0, outdir="./"):
	headers = []
	seqs = []
	i=0
	cum_length=0
	printcount=1
	for line in open(str(fasta)).readlines():
		if line.strip().startswith(">"):
			headers.append(line.strip())
			seqs.append("")
			if i >= 1:
				cum_length+=len(seqs[-2])
#				print("%s\t%s\t%s" %(headers[-2], len(seqs[-2]), cum_length))
			if cum_length >= max_length:
				if pr:
					os.mkdir(outdir+"/"+str(printcount).zfill(4))
					fh = open(outdir+"/"+str(printcount).zfill(4)+"/p0001", 'w')
#					print("%s\t%s" %(str(printcount).zfill(4), cum_length)) #"{:04d}".format(printcount))
					for j in range(len(headers)-1):
						fh.write("%s\n%s\n" %(headers[j],seqs[j]))
					fh.close()
				for j in reversed(range(len(headers)-1)):
					del headers[j]
					del seqs[j]
					cum_length=len(seqs[-1])
#				print("the lenght is again: %s" %len(headers))
				printcount+=1
			i+=1
		else:
			seqs[-1] = seqs[-1]+line.strip()

	if pr:
		os.mkdir(outdir+"/"+str(printcount).zfill(4))
		fh = open(outdir+"/"+str(printcount).zfill(4)+"/p0001", 'w')

#		print("%s\t%s" %(str(printcount).zfill(4), cum_length+len(seqs[-1])))
		for j in range(len(headers)):
			fh.write("%s\n%s\n" %(headers[j],seqs[j]))
		fh.close()

	if not pr:
		return printcount

for sample in samples.index.values.tolist():
    counter=parition_by_length(str(samples.fasta[sample]), max_length=n, pr=0) 
    for i in range(1,counter+1):
        dic['sample'].append(sample)
        dic['unit'].append(str(i).zfill(4))
	

#for sample in samples.index.values.tolist():
##    print sample,samples.fasta[sample]
#    counter=0
#    for line in open(samples.fasta[sample]):
#        if line.startswith(">"):
#            counter+=1
##    print int(ceil(counter/float(n)))
#    for i in range(1,int(ceil(counter/float(n)))+1):
#        dic['sample'].append(sample)
#        dic['unit'].append(str(i).zfill(4))
#
##print dic
units = pd.DataFrame(dic).set_index(['sample','unit'], drop=False)
#print units
#for row in units.itertuples():
#    print(row)

#units = pd.read_csv(config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index

#print(units)
#print(expand("{unit}", unit=units.index.tolist()))

def get_assembly_path(wildcards):
# this is to get the path to the assembly from the CSV file
	return samples.loc[wildcards.sample, ["fasta"]].to_list()

def calculate_n_split(wildcards):
	from math import ceil
	counter=0
	n=1
	for l in open(get_assembly_path(wildcards)).readlines():
		if l.startswith(">"):
			counter+=1
	return list(range(1,ceil(counter/n)+1))


def get_protein_evidence_path(p="data/protein_evidence/*.fasta.gz"):
	import glob, os
	return [os.path.abspath(x) for x in glob.glob(p)]

rule all:
	input:
#		expand("results/{name}/BUSCO/run_{name}/single_copy_busco_sequences/{name}.BUSCOs.fasta", name=samples.index.tolist()),
#		expand("results/{name}/SNAP.PASS1/snap.status.ok", name=samples.index.tolist()),
#		expand("results/{name}/REPEATMASKER/repeatmasker.status.ok", name=samples.index.tolist()),
#		expand("results/{name}/REPEATMASKER/{name}.masked.final.out.reformated.gff", name=samples.index.tolist()),
#		expand("results/{name}/NR_PROTEIN_EVIDENCE/nr.status.ok", name=samples.index.tolist()),
		expand("results/{name}/GENOME_PARTITIONS/splitting.ok", name=samples.index.tolist()),
#		expand("results/{name}/MAKER.PASS1/init.ok", name=samples.index.tolist())
#		expand("results/{name}/MAKER.PASS1/splitting.ok", name=samples.index.tolist())
#		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.fasta", unit=units.itertuples())
#		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.all.maker.gff", unit=units.itertuples()),
#		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.noseq.maker.gff", unit=units.itertuples()),
		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.maker.output.tar.gz", unit=units.itertuples()),
#		expand("results/{name}/MAKER.PASS1/{name}.all.maker.gff", name=samples.index.tolist()),
#		expand("results/{name}/MAKER.PASS1/{name}.noseq.maker.gff", name=samples.index.tolist()),
#		expand("results/{name}/SNAP.PASS2/{name}.MAKER_PASS1.snap.hmm", name=samples.index.tolist()),
#		expand("results/{name}/AUGUSTUS.PASS2/augustus.ok", name=samples.index.tolist())
		expand("results/{name}/MAKER.PASS2/init.ok", name=samples.index.tolist()),
#		expand("results/{name}/MAKER.PASS2/splitting.ok", name=samples.index.tolist()),
		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.all.maker.gff", unit=units.itertuples()),
		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.noseq.maker.gff", unit=units.itertuples()),
		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.maker.output.tar.gz", unit=units.itertuples()),
		expand("results/{name}/MAKER.PASS2/{name}.all.maker.gff", name=samples.index.tolist()),
		expand("results/{name}/MAKER.PASS2/{name}.noseq.maker.gff", name=samples.index.tolist())

rule initiate:
        params:
                prefix = "{sample}"
        output:
                "results/{sample}/{sample}.ok"
        shell:
                """
		if [[ ! -d results/{params.prefix} ]]
		then
			mkdir results/{params.prefix}
		fi
		touch {output}
		"""

rule busco:
	input:
		ok = rules.initiate.output,
		fasta = get_assembly_path
	params:
		prefix = "{sample}",
		busco_set = "data/BUSCO/arthropoda_odb9/",
		augustus_species = "fly"
	threads: 8
	singularity:
		"docker://chrishah/busco-docker:v3.1.0"
	log:
		stdout = "results/{sample}/logs/BUSCO.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/BUSCO.{sample}.stderr.txt"
	output:
		"results/{sample}/BUSCO/run_{sample}/single_copy_busco_sequences/{sample}.BUSCOs.fasta"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/BUSCO ]]
		then
			mkdir results/{params.prefix}/BUSCO
		fi
		cd results/{params.prefix}/BUSCO

		if [[ ! -d tmp ]]
		then
			mkdir tmp
		fi

		cp -rf /usr/share/augustus/config tmp/config
		AUGUSTUS_CONFIG_PATH=$(pwd)/tmp/config

		#run BUSCO
		run_BUSCO.py \
		--in ../../../{input.fasta} --out {params.prefix} -l ../../../{params.busco_set} --mode genome -c {threads} -f \
		-sp {params.augustus_species} --long --augustus_parameters='--progress=true' 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		#collect predicted BUSCOs
		cat run_{params.prefix}/single_copy_busco_sequences/*.faa | sed 's/:.*//' > run_{params.prefix}/single_copy_busco_sequences/{params.prefix}.BUSCOs.fasta

		echo -e "\n$(date)\tFinished!\n"
		"""

rule cegma:
	input:
		ok = rules.initiate.output,
		fasta = get_assembly_path
	params:
		prefix = "{sample}"
	threads: 10
	singularity:
		"docker://chrishah/cegma:2.5"
	log:
		stdout = "results/{sample}/logs/CEGMA.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/CEGMA.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/CEGMA/cegma.status.ok",
		cegma_gff = "results/{sample}/CEGMA/{sample}.cegma.gff"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/CEGMA ]]
		then
			mkdir results/{params.prefix}/CEGMA
		fi
		cd results/{params.prefix}/CEGMA

		#run CEGMA
		cegma -g ../../../{input.fasta} -T {threads} -o {params.prefix} 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		retVal=$?
		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "Cegma ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi

		"""

rule snap_pass1:
	input:
		ok = rules.cegma.output.ok,
		cegma_gff = rules.cegma.output.cegma_gff,
		fasta = get_assembly_path
	params:
		prefix = "{sample}",
		aed = ""
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/SNAP.PASS1.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/SNAP.PASS1.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/SNAP.PASS1/snap.status.ok",
		hmm = "results/{sample}/SNAP.PASS1/{sample}.cegma.snap.hmm"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/SNAP.PASS1 ]]
		then
			mkdir results/{params.prefix}/SNAP.PASS1
		fi
		cd results/{params.prefix}/SNAP.PASS1

		echo -e "[$(date)]\tConvert CEGMA gff to SNAP input"
		cegma2zff ../../../{input.cegma_gff} ../../../{input.fasta}

		echo -e "[$(date)]\tgather some stats and validate"
		fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
		fathom genome.ann genome.dna -validate > validate.log 2>&1
		
		echo -e "[$(date)]\tcollect the training sequences and annotations, plus 1000 surrounding bp for training"
		fathom genome.ann genome.dna -categorize 1000
		fathom -export 1000 -plus uni.ann uni.dna

		echo -e "[$(date)]\tcreate the training parameters"
		forge export.ann export.dna

		echo -e "[$(date)]\tassemble the HMMs"
		hmm-assembler.pl {params.prefix} . > ../../../{output.hmm}

		retVal=$?
		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "SNAP ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		"""
rule repeatmodeler:
	input:
		fasta = get_assembly_path
	params:
		prefix = "{sample}",
	threads: 5
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/REPEATMODELER.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/REPEATMODELER.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/REPEATMODELER/repeatmodeler.status.ok",
		fasta = "results/{sample}/REPEATMODELER/{sample}-families.fa"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/REPEATMODELER ]]
		then
			mkdir results/{params.prefix}/REPEATMODELER
		fi
		cd results/{params.prefix}/REPEATMODELER

		#run REPEATMODELER
		BuildDatabase -name {params.prefix} -engine ncbi ../../../{input.fasta} 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		RepeatModeler -pa {threads} -engine ncbi -database {params.prefix} 1>> ../../../{log.stdout} 2>> ../../../{log.stderr}

		retVal=$?
		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "REPEATMODELER ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		"""


rule repeatmasker:
	input:
		repmod = rules.repeatmodeler.output.fasta,
		fasta = get_assembly_path
	params:
		prefix = "{sample}",
		repeat_taxon = "eukaryota",
		conversion_script = "bin/convert_repeatmasker_gff_to_MAKER_compatible_gff.sh"
	threads: 10
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/REPEATMASKER.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/REPEATMASKER.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/REPEATMASKER/repeatmasker.status.ok",
		gff = "results/{sample}/REPEATMASKER/{sample}.masked.final.out.reformated.gff",
		masked = "results/{sample}/REPEATMASKER/{sample}.masked.final.fasta"
	shell:
		"""
		basedir=$(pwd)

		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/REPEATMASKER ]]
		then
			mkdir results/{params.prefix}/REPEATMASKER
		fi
		cd results/{params.prefix}/REPEATMASKER

		#this is a bit of a hack, but since singularity does not allow to directly write to images ('no space left') and RepeatMasker in the
		#process needs to produce some files, I need to first get the RepeatMasker directory out of the image.
		#Then use the executable in the directory in my writable environment
		#apparently I could also use singularities '--sandbox' option, but from what I can see this would write the content of the entire image to a 
		#directory, so it would take much longer

		#Copy the RepeatMasker directory from the image
		cp -pfr /usr/local/RepeatMasker .

		#Do RepeatMasking with denovo library
		mkdir denovo
		./RepeatMasker/RepeatMasker -engine ncbi -s -pa {threads} -lib $basedir/{input.repmod} -noisy -dir denovo -gff $basedir/{input.fasta} 1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		retVal=$?

		#run REPEATMASKER against full repeat library, but use only the assembly that is already masked based on the denovo library
		mkdir full
		ln -s $(find ./denovo -name '*fasta.masked') {params.prefix}.masked.denovo.fasta
		./RepeatMasker/RepeatMasker -engine ncbi -s -pa {threads} -species {params.repeat_taxon} -noisy -dir full -gff {params.prefix}.masked.denovo.fasta 1>> $basedir/{log.stdout} 2>> $basedir/{log.stderr}
		retVal=$(( retVal + $? ))

		#cleanup - remove the RepeatMasker directory
		rm -rf RepeatMasker
		rm {params.prefix}.masked.denovo.fasta

		#produce the final repeat annotation
		#copy the final masked fasta and out files from the last Repeatmasker run
		mkdir final
		cd final
		ln -s ../full/{params.prefix}.masked.denovo.fasta.masked {params.prefix}.masked.final.fasta
		ln -s ../full/{params.prefix}.masked.denovo.fasta.out {params.prefix}.masked.final.out

		#produce gff3 file from the final RepeatMasker output (this gff3 seems to work well with MAKER after some conversion - see below)
		/usr/local/RepeatMasker/util/rmOutToGFF3.pl {params.prefix}.masked.final.out > {params.prefix}.masked.final.out.gff3
		retVal=$(( retVal + $? ))

		#modify gff3 file so MAKER accepts it down the line
		$basedir/{params.conversion_script} {params.prefix}.masked.final.out.gff3 > $basedir/{output.gff}

		cd ..
		ln -s final/{params.prefix}.masked.final.fasta $basedir/{output.masked}
		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error"
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		"""

rule prepare_protein_evidence:
	input:
		expand("{full}/{file}", full=[os.getcwd()], file=glob.glob("data/protein_evidence/*.fasta.gz"))
#		get_protein_evidence_path
#		["/home/lv71312/hahnc/SNAKEMAKE/maker/data/protein_evidence/uniprot_sprot.fasta.gz", "/home/lv71312/hahnc/SNAKEMAKE/maker/data/protein_evidence/other.fasta.gz"]
	params:
		prefix = "{sample}",
		mem = "8000",
		similarity = "0.98"
	threads: 2
	singularity:
		"docker://chrishah/cdhit:v4.8.1"
	log:
		stdout = "results/{sample}/logs/CDHIT.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/CDHIT.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/NR_PROTEIN_EVIDENCE/nr.status.ok",
		nr_proteins = "results/{sample}/NR_PROTEIN_EVIDENCE/nr_external_proteins.cd-hit.fasta"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/NR_PROTEIN_EVIDENCE ]]
		then
			mkdir results/{params.prefix}/NR_PROTEIN_EVIDENCE
		fi
		cd results/{params.prefix}/NR_PROTEIN_EVIDENCE

		echo -e "Remove redundancy at {params.similarity} in files:\n{input}"

		#concatenate all physical evidence
		cat {input} > external_proteins.fasta.gz

		#run cd-hit
		cd-hit -T {threads} -M {params.mem} -i external_proteins.fasta.gz -o external_proteins.cd-hit-{params.similarity}.fasta -c {params.similarity} 1> ../../../{log.stdout} 2> ../../../{log.stderr}
		retVal=$?

		ln -s external_proteins.cd-hit-{params.similarity}.fasta ../../../{output.nr_proteins}

		#remove obsolete file
		rm external_proteins.fasta.gz

		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi

		"""

rule split_fasta:
	input:
		fasta = get_assembly_path
	params:
		prefix = "{sample}",
		outdir = "results/{sample}/GENOME_PARTITIONS/"
#		len = n
	output:
		ok = "results/{sample}/GENOME_PARTITIONS/splitting.ok"
	run:
		parition_by_length(input.fasta, max_length=n, pr=1, outdir=params.outdir)
		Path(output.ok).touch()		

rule initiate_MAKER_PASS1:
	input:
		snap = rules.snap_pass1.output.hmm,
		nr_evidence = rules.prepare_protein_evidence.output.nr_proteins,
		busco_proteins = rules.busco.output,
		repmod_lib = rules.repeatmodeler.output.fasta,
		repmas_gff = rules.repeatmasker.output.gff
	params:
		prefix = "{sample}",
		alt_est = expand("{full}/{files}", full=[os.getcwd()], files=glob.glob("data/transcripts_alt/*")),
		est = expand("{full}/{files}", full=[os.getcwd()], files=glob.glob("data/transcripts/*"))
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS1.init.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS1.init.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/MAKER.PASS1/init.ok"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/MAKER.PASS1 ]]
		then
			mkdir results/{params.prefix}/MAKER.PASS1
		fi
		cd results/{params.prefix}/MAKER.PASS1

		maker -CTL 1> ../../../{log.stdout} 2> ../../../{log.stderr}
		retVal=$?

		##### Modify maker_opts.ctl file
		#add SNAP result if present
		if [ -f "../../../{input.snap}" ]
		then
        		echo -e "SNAP hmms provided: {input.snap}"
        		sed -i "s?snaphmm= ?snaphmm=$(pwd)/../../../{input.snap} ?" maker_opts.ctl
		fi
		#add protein evidence if present
		if [ -f "../../../{input.nr_evidence}" ] || [ -f "../../../{input.busco_proteins}" ]
		then
			paths=""
        		if [ -f "../../../{input.nr_evidence}" ]; then echo -e "Protein evidence provided: {input.nr_evidence}"; paths="$paths,$(pwd)/../../../{input.nr_evidence}"; fi
        		if [ -f "../../../{input.busco_proteins}" ]; then echo -e "BUSCO proteins provided: {input.busco_proteins}"; paths="$paths,$(pwd)/../../../{input.nr_evidence}"; fi
			paths=$(echo $paths | sed 's/^,//')
        		sed -i "s?protein= ?protein=$paths ?" maker_opts.ctl
		fi
		#add denovo repeat library if present
		if [ -f "../../../{input.repmod_lib}" ]
		then
		        echo -e "Denovo repeat library provided: {input.repmod_lib}"
		        sed -i "s?rmlib= ?rmlib=$(pwd)/../../../{input.repmod_lib} ?" maker_opts.ctl
		fi
		#add repeatmasker gff if present
		if [ -f "../../../{input.repmas_gff}" ]
		then
        		echo -e "Repeatmasker gff provided: {input.repmas_gff}"
        		sed -i "s?rm_gff= ?rm_gff=$(pwd)/../../../{input.repmas_gff} ?" maker_opts.ctl
		fi

		alt_est="{params.alt_est}"
		if [ ! -z "$alt_est" ]
		then
			echo -e "Transcriptome evidence (altest) provided: {params.alt_est}"
			sed -i "s?altest= ?altest=$(echo $alt_est | sed 's/ /,/g') ?" maker_opts.ctl
		fi

		est="{params.est}"
		if [ ! -z "$est" ]
		then
			echo -e "Transcriptome evidence (est) provided: {params.est}"
			sed -i "s?^est= ?est=$(echo $est | sed 's/ /,/g') ?" maker_opts.ctl
		fi

		sed -i 's/est2genome=0/est2genome=1/' maker_opts.ctl
		sed -i 's/protein2genome=0/protein2genome=1/' maker_opts.ctl
		############################

		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		"""

#rule split_fasta_PASS1:
#	input:
#		rules.initiate_maker_pass1.output.ok,
#		fasta = get_assembly_path
#	params:
#		prefix = "{sample}",
#		seqs_per_file = "1"
#	singularity:
#		"docker://chrishah/ectools-docker"
#	output:
#		ok = "results/{sample}/MAKER.PASS1/splitting.ok"
##		file = "results/{sample}/MAKER.PASS1/{unit}/p0001"
#	shell:
#		"""
#		#get going
#		echo -e "\n$(date)\tStarting ...\n"
#
#		cd results/{params.prefix}/MAKER.PASS1
#		
#		echo -e "\nSplitting up assembly ({input.fasta}): {params.seqs_per_file} sequence(s) per file\n"
#		partition.py {params.seqs_per_file} 1 ../../../{input.fasta}
#		
#		retVal=$?
#
#		echo -e "\n$(date)\tFinished!\n"
#
#		if [ ! $retVal -eq 0 ]
#		then
#			echo "There was some error"
#			exit $retVal
#		else
#			touch ../../../{output.ok}
#		fi
#		"""

rule run_MAKER_PASS1:
	input:
		init_ok = rules.initiate_MAKER_PASS1.output.ok,
		split_ok = rules.split_fasta.output.ok
	params:
		sub = "results/{sample}/GENOME_PARTITIONS/{unit}/p0001",
		dir = "{unit}",
		prefix = "{sample}"
	threads: 9
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS1.run.{sample}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS1.run.{sample}.{unit}.stderr.txt"
	output:
		sub_fasta = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.fasta",
		gff = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.noseq.maker.gff"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		basedir=$(pwd)

		cd results/{params.prefix}/MAKER.PASS1/{params.dir}
		ln -s $basedir/{params.sub} {params.prefix}.{params.dir}.fasta

		#run MAKER
		maker -base {params.prefix}.{params.dir} -g {params.prefix}.{params.dir}.fasta -nolock -c {threads} ../maker_opts.ctl ../maker_bopts.ctl ../maker_exe.ctl 1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		#prepare data from MAKER 
		cd {params.prefix}.{params.dir}.maker.output
		gff3_merge -d {params.prefix}.{params.dir}_master_datastore_index.log -o {params.prefix}.{params.dir}.all.maker.gff
		fasta_merge -d {params.prefix}.{params.dir}_master_datastore_index.log
		gff3_merge -n -d {params.prefix}.{params.dir}_master_datastore_index.log -o {params.prefix}.{params.dir}.noseq.maker.gff
		cd ..

		mv {params.prefix}.{params.dir}.maker.output/{params.prefix}.{params.dir}.all.maker.* .
		mv {params.prefix}.{params.dir}.maker.output/{params.prefix}.{params.dir}.noseq.maker.* .
		
		echo -e "\n$(date)\tFinished!\n"
		"""

rule cleanup_MAKER_PASS1:
	input:
		rules.run_MAKER_PASS1.output
	params:
		dir = "{unit}",
		prefix = "{sample}",
		script = "bin/cleanup.sh"
	output:
		"results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.maker.output.tar.gz"
	shell:
		"""
		basedir=$(pwd)
		
		cd results/{params.prefix}/MAKER.PASS1/{params.dir}/
		bash $basedir/{params.script} {params.prefix}.{params.dir}.maker.output

		"""	

rule merge_MAKER_PASS1:
	input:
		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.all.maker.gff", unit=units.itertuples())
	params:
		prefix = "{sample}",
		script = "bin/merging.sh"
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	output:
		all_gff = "results/{sample}/MAKER.PASS1/{sample}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.gff",
		proteins = "results/{sample}/MAKER.PASS1/{sample}.all.maker.proteins.fasta"
	shell:
		"""
		basedir=$(pwd)
		cd results/{params.prefix}/MAKER.PASS1/

		bash $basedir/{params.script} {params.prefix}
		"""

rule snap_pass2:
	input:
		rules.merge_MAKER_PASS1.output.all_gff
	params:
		aed = "config["aed"]["snap_pass2"],
		prefix = "{sample}",
		script = "bin/snap.p2.sh"
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/SNAP.PASS2.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/SNAP.PASS2.{sample}.stderr.txt"
	output:
		snap_hmm = "results/{sample}/SNAP.PASS2/{sample}.MAKER_PASS1.snap.hmm"
	shell:
		"""
		basedir=$(pwd)
		if [[ ! -d results/{params.prefix}/SNAP.PASS2 ]]
		then
			mkdir results/{params.prefix}/SNAP.PASS2/
		fi
		cd results/{params.prefix}/SNAP.PASS2/

		bash $basedir/{params.script} \
		{params.prefix} \
		$basedir/{input} \
		{params.aed} \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		"""

rule AUGUSTUS_PASS2:
	input:
		busco_ok = rules.busco.output,
		fasta = rules.repeatmasker.output.masked,
		maker_proteins = rules.merge_MAKER_PASS1.output.proteins
	params:
		prefix = "{sample}",
		training_params = "results/{sample}/BUSCO/run_{sample}/augustus_output/retraining_parameters",
		script = "bin/augustus.PASS2.sh",
		aed = "config["aed"]["AUGUSTUS_PASS2"]
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/AUGUSTUS.PASS2.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/AUGUSTUS.PASS2.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/AUGUSTUS.PASS2/augustus.ok",
		abinitio = "results/{sample}/AUGUSTUS.PASS2/augustus.gff3",
		training_params = directory("results/{sample}/AUGUSTUS.PASS2/{sample}")
	shell:
		"""
		basedir=$(pwd)
		
		echo -e "TODO: CHECK FOR cDNA evidence and include in autoAug.pl run via --cdna=cdna.fa option - see 'Typical Usage' in Readme of autoAug.pl script"		
		echo -e "TODO: RUN Augustus across a range of aed cutoffs and use the one that has the best prediction accuracy"		

		if [[ ! -d results/{params.prefix}/AUGUSTUS.PASS2 ]]
		then
			mkdir results/{params.prefix}/AUGUSTUS.PASS2/
		fi
		cd results/{params.prefix}/AUGUSTUS.PASS2/

		if [[ ! -d tmp ]]
		then
			mkdir tmp
		fi
		
		#copy augustus config directory from the image
		cp -rf /usr/share/augustus/config tmp/config

		bash $basedir/{params.script} \
		{params.prefix} \
		$basedir/{input.fasta} \
		$basedir/{input.maker_proteins} \
		{params.aed} \
		$(pwd)/tmp/config \
		$basedir/{params.training_params} \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		retVal=$?

		rm -rf tmp/

		if [ ! $retVal -eq 0 ]
		then
			>&2 echo "Augustus ended in an error" >> $basedir/{log.stderr}
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi

		"""

rule initiate_MAKER_PASS2:
	input:
		pred_gff = rules.AUGUSTUS_PASS2.output.abinitio,
		params = rules.AUGUSTUS_PASS2.output.training_params,
		snaphmm = rules.snap_pass2.output.snap_hmm,
		MP1_ok = rules.merge_MAKER_PASS1.output
	params:
		prefix = "{sample}",
		script = "bin/prepare_maker_opts_PASS2.sh",
		protein_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.protein2genome.gff",
		rm_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.repeats.gff",
		altest_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.cdna2genome.gff",
		est_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.est2genome.gff"
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS2.init.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS2.init.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/MAKER.PASS2/init.ok"
	shell:
		"""
		basedir=$(pwd)

		#get going
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/MAKER.PASS2 ]]
		then
			mkdir results/{params.prefix}/MAKER.PASS2
		fi
		cd results/{params.prefix}/MAKER.PASS2

		maker -CTL 1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		retVal=$?

		if [[ ! -d tmp ]]
		then
			mkdir tmp
		fi
		#copy augustus config directory from the image
		cp -rf /usr/share/augustus/config tmp/config

		##### Modify maker_opts.ctl file
		bash $basedir/{params.script} \
		$basedir/{input.snaphmm} \
		{params.prefix} \
		$basedir/{input.params} \
		$basedir/{input.pred_gff} \
		$basedir/{params.rm_gff} \
		$basedir/{params.protein_gff} \
		$basedir/{params.altest_gff} \
		$basedir/{params.est_gff} \
		$(pwd)/tmp/config \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		
		retVal=$(( retVal + $? ))

		echo -e "\n$(date)\tFinished!\n"

		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error" >> $basedir/{log.stderr}
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		
		"""

#rule split_fasta_PASS2:
#	input:
#		rules.initiate_MAKER_PASS2.output.ok,
#		fasta = get_assembly_path
#	params:
#		prefix = "{sample}",
#		seqs_per_file = "1"
#	singularity:
#		"docker://chrishah/ectools-docker"
#	output:
#		ok = "results/{sample}/MAKER.PASS2/splitting.ok"
##		file = "results/{sample}/MAKER.PASS2/{unit}/p0001"
#	shell:
#		"""
#		basedir=$(pwd)
#
#		#get going
#		echo -e "\n$(date)\tStarting ...\n"
#
#		cd results/{params.prefix}/MAKER.PASS2
#		
#		echo -e "\nSplitting up assembly ({input.fasta}): {params.seqs_per_file} sequence(s) per file\n"
#		partition.py {params.seqs_per_file} 1 $basedir/{input.fasta}
#		
#		retVal=$?
#
#		echo -e "\n$(date)\tFinished!\n"
#
#		if [ ! $retVal -eq 0 ]
#		then
#			echo "There was some error"
#			exit $retVal
#		else
#			touch $basedir/{output.ok}
#		fi
#		"""


rule run_MAKER_PASS2:
	input:
		split_ok = rules.split_fasta.output.ok,
		init_ok = rules.initiate_MAKER_PASS2.output.ok
	params:
		sub = "results/{sample}/GENOME_PARTITIONS/{unit}/p0001",
		dir = "{unit}",
		prefix = "{sample}"
	threads: 2
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS2.run.{sample}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS2.run.{sample}.{unit}.stderr.txt"
	output:
		sub_fasta = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.fasta",
		gff = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.noseq.maker.gff"
	shell:
		"""
		#get going
		echo -e "\n$(date)\tStarting ...\n"

		basedir=$(pwd)

		cd results/{params.prefix}/MAKER.PASS2/{params.dir}
		ln -s $basedir/{params.sub} {params.prefix}.{params.dir}.fasta

		AUGUSTUS_CONFIG_PATH=$basedir/results/{params.prefix}/MAKER.PASS2/tmp/config

		#run MAKER
		maker -base {params.prefix}.{params.dir} -g {params.prefix}.{params.dir}.fasta -nolock -c {threads} ../maker_opts.ctl ../maker_bopts.ctl ../maker_exe.ctl 1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		#prepare data from MAKER 
		cd {params.prefix}.{params.dir}.maker.output
		gff3_merge -d {params.prefix}.{params.dir}_master_datastore_index.log -o {params.prefix}.{params.dir}.all.maker.gff
		fasta_merge -d {params.prefix}.{params.dir}_master_datastore_index.log
		gff3_merge -n -d {params.prefix}.{params.dir}_master_datastore_index.log -o {params.prefix}.{params.dir}.noseq.maker.gff
		cd ..

		mv {params.prefix}.{params.dir}.maker.output/{params.prefix}.{params.dir}.all.maker.* .
		mv {params.prefix}.{params.dir}.maker.output/{params.prefix}.{params.dir}.noseq.maker.* .
		
		echo -e "\n$(date)\tFinished!\n"
		"""

rule cleanup_MAKER_PASS2:
	input:
		rules.run_MAKER_PASS2.output
	params:
		dir = "{unit}",
		prefix = "{sample}",
		script = "bin/cleanup.sh"
	output:
		"results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.maker.output.tar.gz"
	shell:
		"""
		basedir=$(pwd)
		
		cd results/{params.prefix}/MAKER.PASS2/{params.dir}/

		bash $basedir/{params.script} {params.prefix}.{params.dir}.maker.output

		"""	

rule merge_MAKER_PASS2:
	input:
		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.all.maker.gff", unit=units.itertuples())
	params:
		prefix = "{sample}",
		script = "bin/merging.sh"
	singularity:
		"docker://chrishah/maker-full:2.31.10"
	output:
		all_gff = "results/{sample}/MAKER.PASS2/{sample}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS2/{sample}.noseq.maker.gff",
		proteins = "results/{sample}/MAKER.PASS2/{sample}.all.maker.proteins.fasta"
	shell:
		"""
		basedir=$(pwd)
		cd results/{params.prefix}/MAKER.PASS2/

		bash $basedir/{params.script} {params.prefix}
		"""
