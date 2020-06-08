configfile: "data/config.yaml"

import pandas as pd
import os
import glob
from math import ceil
from pathlib import Path
from subprocess import call

localrules: initiate, all

n=int(config["split_batch_length"])
min=int(config["split_min_length"])

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]

dic = {'sample': [], 'unit': []}

def get_assembly_path(wildcards):
# this is to get the path to the assembly from the CSV file
	return samples.loc[wildcards.sample, ["fasta"]].to_list()

def get_transcripts_path(wildcards, p="data/transcripts/*"):
	#get paths to fasta transcript fasta files - if file has prefix identical to sample prefix in data.csv -> assume it's a transcriptome of this species -> MAKER 'est' option
	dic = {'alt_ests': [], 'ests': []}
	for f in glob.glob(p):
		if f.split("/")[-1].startswith(wildcards.sample):
			dic['ests'].append(os.path.abspath(f))
		else:
			dic['alt_ests'].append(os.path.abspath(f))
	return dic

def partition_by_length(fasta, max_length=n, min_length=min, pr=0, outdir="./"):
#function that partitions the fasta file
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
				if len(seqs[-2]) >= min_length:
					cum_length+=len(seqs[-2])
#					print("%s\t%s\t%s" %(headers[-2], len(seqs[-2]), cum_length))
				else:
					del headers[-2]
					del seqs[-2]
			if cum_length >= max_length:
				if pr:
					if not os.path.exists(outdir+"/"+str(printcount).zfill(4)):
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
		if not os.path.exists(outdir+"/"+str(printcount).zfill(4)):
			os.mkdir(outdir+"/"+str(printcount).zfill(4))
		fh = open(outdir+"/"+str(printcount).zfill(4)+"/p0001", 'w')

#		print("%s\t%s" %(str(printcount).zfill(4), cum_length+len(seqs[-1])))
		for j in range(len(headers)):
			fh.write("%s\n%s\n" %(headers[j],seqs[j]))
		fh.close()

	if not pr:
		return printcount

unitdict = {}
print("Counting partitions (batchsize >= "+str(n)+"bp, minimum length = "+str(min)+"bp) ..")
for sample in samples.index.values.tolist():
    print("\t"+sample+" - n=", end='')
    count = subprocess.run("bash ./bin/count_length.sh %s %i %i count" %(samples.fasta[sample], n, min), shell=True, stdout=subprocess.PIPE)
    counter = int(count.stdout.decode('utf-8').split("\t")[-1])


#    counter=partition_by_length(str(samples.fasta[sample]), max_length=n, min_length=min, pr=0) 
    print(counter)
    print("\t"+count.stdout.decode('utf-8').split("\t")[0])
    unitdict[sample] = []
    for i in range(1,counter+1):
        dic['sample'].append(sample)
        dic['unit'].append(str(i).zfill(4))
	unitdict[sample].append(str(i).zfill(4))	
	
#print(unitdict)
##print dic

units = pd.DataFrame(dic).set_index(['sample','unit'], drop=False)
#print(units)
#print(units.index.tolist())
#print units
#for row in units.itertuples():
#    print(row)

units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index

rule all:
	input:
		expand("results/{unit.sample}/MAKER.PASS1/{unit.unit}/{unit.sample}.{unit.unit}.maker.output.tar.gz", unit=units.itertuples()),
		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.maker.output.tar.gz", unit=units.itertuples()),
		expand("results/{name}/REPEATMODELER/repeatmodeler.cleanup.ok", name=samples.index.tolist()),
		expand("results/{name}/MAKER.PASS2/{name}.all.maker.gff", name=samples.index.tolist())

rule setup_maker:
	input:
		maker_tarball = config["maker_tarball"],
	params:
		repbase = config["RepbaseRepeatMaskerEdition"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	output: 
		bin = directory("bin/maker/bin")
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)

		#Copy the RepeatMasker directory from the image
		cp -pfrv /usr/local/RepeatMasker bin/

		if [ "{params.repbase}" == "None" ]
		then
			echo -e "No additional Repeatlibrary provided - ok"
		else
			bin/setup_Repeatmasker.sh bin/ {params.repbase}
#			tar xvfz {params.repbase} -C bin/RepeatMasker/
#			perl bin/RepeatMasker/rebuild
		fi

		if [ "{input.maker_tarball}" == "None" ]
		then
			echo -e "Providing a maker tarball is mandatory"
			exit 1
		else
			bash bin/setup_maker.sh {input.maker_tarball} bin 
		fi


		echo -e "\n$(date)\tFinished!\n"
		"""
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

rule split:
	input:
		fasta = get_assembly_path,
		ok = rules.initiate.output
	params:
		prefix = "{sample}",
		len = n,
		min = min
	log:
		stdout = "results/{sample}/logs/split.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/split.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/GENOME_PARTITIONS/splitting.ok",
		fasta = "results/{sample}/GENOME_PARTITIONS/{sample}.min"+str(min)+".fasta"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)

		cd results/{params.prefix}/GENOME_PARTITIONS/
		bash $basedir/bin/count_length.sh ../../../{input.fasta} {params.len} {params.min} split

		retVal=$?

                if [ ! $retVal -eq 0 ]
                then
                        echo "Splitting ended in an error"
                        exit $retVal
                else
                        touch ../../../{output.ok}
			cat *.fasta > ../../../{output.fasta}
                fi

		echo -e "\n$(date)\tFinished!\n"
		"""

rule genemark:
	input:
		ok = rules.initiate.output,
		fasta = rules.split.output.fasta
	params:
		prefix = "{sample}",
		genemark_dir = config["genemark"]["genemark_dir"],
		gmes_petap_params = config["genemark"]["gmes_petap_params"]
	threads: config["threads"]["genemark"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/GENEMARK.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/GENEMARK.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/GENEMARK/genemark.status.ok",
		model = "results/{sample}/GENEMARK/gmhmm.mod"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)

                if [[ ! -d results/{params.prefix}/GENEMARK ]]
                then
                        mkdir results/{params.prefix}/GENEMARK
		else
			if [ "$(ls -1 results/{params.prefix}/GENEMARK/ | wc -l)" -gt 0 ]
			then
				echo -e "Cleaning up remnants of previous run first" 1> {log.stdout} 2> {log.stderr}
				rm results/{params.prefix}/GENEMARK
				mkdir results/{params.prefix}/GENEMARK
			fi
                fi
                cd results/{params.prefix}/GENEMARK

		ln -sf $basedir/{params.genemark_dir}/gm_key .gm_key

		if [ "{params.gmes_petap_params}" == "None" ]
		then
			gmes_petap.pl -ES -cores {threads} -sequence ../../../{input.fasta} 1> ../../../{log.stdout} 2> ../../../{log.stderr}
		else
			gmes_petap.pl -ES {params.gmes_petap_params} -cores {threads} -sequence ../../../{input.fasta} 1> ../../../{log.stdout} 2> ../../../{log.stderr}
		fi

		retVal=$?

		if [ ! $retVal -eq 0 ]
		then
			echo "Genemark ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		
		"""		
		
rule busco:
	input:
		ok = rules.initiate.output,
		fasta = rules.split.output.fasta
	params:
		prefix = "{sample}",
		busco_path = "data/BUSCO",
		busco_set = config["busco"]["set"],
		augustus_species = config["busco"]["species"]
	threads: config["threads"]["busco"]
	singularity:
		"docker://chrishah/busco-docker:v3.1.0"
	log:
		stdout = "results/{sample}/logs/BUSCO.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/BUSCO.{sample}.stderr.txt"
	output:
		"results/{sample}/BUSCO/run_{sample}/single_copy_busco_sequences/{sample}.BUSCOs.fasta"
	shell:
		"""
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
		--in ../../../{input.fasta} --out {params.prefix} -l ../../../{params.busco_path}/{params.busco_set} --mode genome -c {threads} -f \
		-sp {params.augustus_species} --long --augustus_parameters='--progress=true' 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		#collect predicted BUSCOs
		cat run_{params.prefix}/single_copy_busco_sequences/*.faa | sed 's/:.*//' > run_{params.prefix}/single_copy_busco_sequences/{params.prefix}.BUSCOs.fasta

		echo -e "\n$(date)\tFinished!\n"
		"""

rule cegma:
	input:
		ok = rules.initiate.output,
		fasta = rules.split.output.fasta
	params:
		prefix = "{sample}"
	threads: config["threads"]["cegma"]
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
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/CEGMA ]]
		then
			mkdir results/{params.prefix}/CEGMA
		fi
		cd results/{params.prefix}/CEGMA

		#run CEGMA
		cegma -g ../../../{input.fasta} -T {threads} -o {params.prefix} 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		retVal=$?

		if [ ! $retVal -eq 0 ]
		then
			echo "Cegma ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"

		"""

rule snap_pass1:
	input:
		ok = rules.cegma.output.ok,
		cegma_gff = rules.cegma.output.cegma_gff,
		fasta = rules.split.output.fasta
	params:
		prefix = "{sample}",
		script = "bin/snap.p1.sh"
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/SNAP.PASS1.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/SNAP.PASS1.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/SNAP.PASS1/snap.status.ok",
		hmm = "results/{sample}/SNAP.PASS1/{sample}.cegma.snap.hmm"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		export PATH="$(pwd)/bin/maker/bin:$PATH"

		if [[ ! -d results/{params.prefix}/SNAP.PASS1 ]]
		then
			mkdir results/{params.prefix}/SNAP.PASS1
		fi
		cd results/{params.prefix}/SNAP.PASS1

		bash $basedir/{params.script} \
		{params.prefix} \
		$basedir/{input.cegma_gff} \
		$basedir/{input.fasta} \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		retVal=$?

		if [ ! $retVal -eq 0 ]
		then
			echo "SNAP ended in an error"
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		"""

rule repeatmodeler:
	input:
		ok = rules.initiate.output,
		fasta = rules.split.output.fasta
	params:
		prefix = "{sample}",
	threads: config["threads"]["repeatmodeler"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/REPEATMODELER.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/REPEATMODELER.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/REPEATMODELER/repeatmodeler.status.ok",
		fasta = "results/{sample}/REPEATMODELER/{sample}-families.fa"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		if [[ ! -d results/{params.prefix}/REPEATMODELER ]]
		then
			mkdir results/{params.prefix}/REPEATMODELER
		else
			if [ "$(ls -1 results/{params.prefix}/REPEATMODELER/ | wc -l)" -gt 0 ]
			then
				echo -e "Cleaning up remnants of previous run first"
				rm results/{params.prefix}/REPEATMODELER
				mkdir results/{params.prefix}/REPEATMODELER
			fi
		fi
		cd results/{params.prefix}/REPEATMODELER

		#run REPEATMODELER
		BuildDatabase -name {params.prefix} -engine ncbi ../../../{input.fasta} 1> ../../../{log.stdout} 2> ../../../{log.stderr}

		RepeatModeler -pa {threads} -engine ncbi -database {params.prefix} 1>> ../../../{log.stdout} 2>> ../../../{log.stderr}

		retVal=$?

		if [ ! $retVal -eq 0 ]
		then
			echo "REPEATMODELER ended in an error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		"""

rule cleanup_repeatmodeler:
	input:
		rules.repeatmodeler.output
	params:
		prefix = "{sample}"
	log:
		stdout = "results/{sample}/logs/REPEATMODELER.cleanup.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/REPEATMODELER.cleanup.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/REPEATMODELER/repeatmodeler.cleanup.ok"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		cd results/{params.prefix}/REPEATMODELER
		for f in $(find ./ -type d -name "RM_*")
		do
			echo -e "\nCompressing $f\n"
			tar cfz $f.tar.gz $f
			
			if [ $? -eq 0 ]
			then
        			rm -rf $f
			else
        			echo -e "Some problem with $f"
			fi
		done

		cd $basedir
		touch {output.ok}

		echo -e "\n$(date)\tFinished!\n"

		"""


rule repeatmasker:
	input:
		fasta = rules.split.output.fasta,
		repmod = rules.repeatmodeler.output.fasta
	params:
		prefix = "{sample}",
		repeat_taxon = "eukaryota",
		conversion_script = "bin/convert_repeatmasker_gff_to_MAKER_compatible_gff.sh"
	threads: config["threads"]["repeatmasker"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/REPEATMASKER.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/REPEATMASKER.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/REPEATMASKER/repeatmasker.status.ok",
		gff = "results/{sample}/REPEATMASKER/{sample}.masked.final.out.reformated.gff",
		masked = "results/{sample}/REPEATMASKER/{sample}.masked.final.fasta"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)

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

		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error"
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		"""

rule prepare_protein_evidence:
	input:
		proteins = expand("{full}/{file}", full=[os.getcwd()], file=glob.glob("data/protein_evidence/*.gz")),
		ok = rules.initiate.output
	params:
		prefix = "{sample}",
		mem = "8000",
		similarity = config["cdhit"]["similarity"]
	threads: config["threads"]["prepare_protein_evidence"]
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

		echo -e "Remove redundancy at {params.similarity} in files:\n{input.proteins}"

		#concatenate all physical evidence
		cat {input.proteins} > external_proteins.fasta.gz

		#run cd-hit
		cd-hit -T {threads} -M {params.mem} -i external_proteins.fasta.gz -o external_proteins.cd-hit-{params.similarity}.fasta -c {params.similarity} 1> ../../../{log.stdout} 2> ../../../{log.stderr}
		retVal=$?

		ln -s external_proteins.cd-hit-{params.similarity}.fasta ../../../{output.nr_proteins}

		#remove obsolete file
		rm external_proteins.fasta.gz


		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error"
			exit $retVal
		else
			touch ../../../{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		"""

rule initiate_MAKER_PASS1:
	input:
		maker = rules.setup_maker.output.bin,
		ok = rules.initiate.output,
		snap = rules.snap_pass1.output.hmm,
		nr_evidence = rules.prepare_protein_evidence.output.nr_proteins,
		busco_proteins = rules.busco.output,
		repmod_lib = rules.repeatmodeler.output.fasta,
		repmas_gff = rules.repeatmasker.output.gff
	params:
		prefix = "{sample}",
		transcripts = get_transcripts_path,
		script = "bin/prepare_maker_opts_PASS1.sh"
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS1.init.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS1.init.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/MAKER.PASS1/init.ok"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		export PATH="$(pwd)/bin/maker/bin:$PATH"


		if [[ ! -d results/{params.prefix}/MAKER.PASS1 ]]
		then
			mkdir results/{params.prefix}/MAKER.PASS1
		fi
		cd results/{params.prefix}/MAKER.PASS1

		maker -CTL 1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		retVal=$?

		##### Modify maker_opts.ctl file
		bash $basedir/{params.script} \
		$basedir/{input.snap} \
		$basedir/{input.nr_evidence} \
		$basedir/{input.busco_proteins} \
		$basedir/{input.repmod_lib} \
		$basedir/{input.repmas_gff} \
		altest="{params.transcripts[alt_ests]}" \
		est="{params.transcripts[ests]}" \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		
		retVal=$(( retVal + $? ))


		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error" >> $basedir/{log.stderr}
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"
		"""

rule run_MAKER_PASS1:
	input:
		maker = rules.setup_maker.output.bin,
		init_ok = rules.initiate_MAKER_PASS1.output.ok,
		split_ok = rules.split.output.ok
	params:
		sub = "results/{sample}/GENOME_PARTITIONS/{unit}.fasta",
		dir = "{unit}",
		prefix = "{sample}"
	threads: config["threads"]["run_MAKER_PASS1"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS1.run.{sample}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS1.run.{sample}.{unit}.stderr.txt"
	output:
		sub_fasta = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.fasta",
		gff = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS1/{unit}/{sample}.{unit}.noseq.maker.gff"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		basedir=$(pwd)
		export PATH="$(pwd)/bin/maker/bin:$PATH"

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
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		cd results/{params.prefix}/MAKER.PASS1/{params.dir}/
		bash $basedir/{params.script} {params.prefix}.{params.dir}.maker.output
		echo -e "\n$(date)\tFinished!\n"

		"""	

rule merge_MAKER_PASS1:
	input:
		lambda wildcards: expand("results/{{sample}}/MAKER.PASS1/{unit}/{{sample}}.{unit}.all.maker.gff", sample=wildcards.sample, unit=unitdict[wildcards.sample])
#		expand("results/{{sample}}/MAKER.PASS1/{unit.unit}/{{sample}}.{unit.unit}.fasta", sample=samples.index.tolist(), unit=units.itertuples())
	params:
		prefix = "{sample}",
		script = "bin/merging.sh"
	singularity:
		"docker://chrishah/premaker-plus:18"
	output:
		all_gff = "results/{sample}/MAKER.PASS1/{sample}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.gff",
		proteins = "results/{sample}/MAKER.PASS1/{sample}.all.maker.proteins.fasta"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		export PATH="$(pwd)/bin/maker/bin:$PATH"
		basedir=$(pwd)
		cd results/{params.prefix}/MAKER.PASS1/

		bash $basedir/{params.script} {params.prefix}
		echo -e "\n$(date)\tFinished!\n"
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
		aed = "{aed}",
		transcripts = get_transcripts_path 
	threads: config["threads"]["AUGUSTUS_PASS2"]
	singularity:
		"docker://chrishah/augustus:v3.3.2"
	log:
		stdout = "results/{sample}/logs/AUGUSTUS.PASS2.{aed}.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/AUGUSTUS.PASS2.{aed}.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/AUGUSTUS.PASS2/{aed}/{aed}.augustus.done",
		training_params = directory("results/{sample}/AUGUSTUS.PASS2/{aed}/{aed}.{sample}")
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		echo -e "TODO: CHECK FOR cDNA evidence and include in autoAug.pl run via --cdna=cdna.fa option - see 'Typical Usage' in Readme of autoAug.pl script"		
		echo -e "TODO: RUN Augustus across a range of aed cutoffs and use the one that has the best prediction accuracy"		

		if [[ ! -d results/{params.prefix}/AUGUSTUS.PASS2 ]]
		then
			mkdir results/{params.prefix}/AUGUSTUS.PASS2/
		fi
		cd results/{params.prefix}/AUGUSTUS.PASS2/

		if [[ ! -d {params.aed} ]]
		then
			mkdir {params.aed}
		fi
		cd {params.aed}

		if [[ ! -d tmp.{params.aed} ]]
		then
			mkdir tmp.{params.aed}
		fi
		
		#copy augustus config directory from the image
		cp -rf /usr/share/augustus/config tmp.{params.aed}/config

		est="{params.transcripts[ests]}"
		if [ ! -z "$est" ]
		then
			cat $est > cdna.{params.aed}.fasta	
		fi

		bash $basedir/{params.script} \
		{threads} \
		{params.aed}.{params.prefix} \
		$basedir/{input.fasta} \
		$basedir/{input.maker_proteins} \
		{params.aed} \
		$(pwd)/tmp.{params.aed}/config \
		$basedir/{params.training_params} \
		cdna.{params.aed}.fasta \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		retVal=$?

		#rm -rf tmp.{params.aed}/

		if [ ! $retVal -eq 0 ]
		then
			>&2 echo "Augustus ended in an error" >> $basedir/{log.stderr}
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		echo -e "\n$(date)\tFinished!\n"

		"""
rule pick_augustus_training_set:
	input:
		lambda wildcards: expand("results/{{sample}}/AUGUSTUS.PASS2/{aed}/{aed}.augustus.done", sample=wildcards.sample, aed=config["aed"]["AUGUSTUS_PASS2"])
	params:
		aeds = expand("{aed}", aed=config["aed"]["AUGUSTUS_PASS2"]),
		prefix = "{sample}"
	output:
		best_params = directory("results/{sample}/AUGUSTUS.PASS2/training_params"),
		gff = "results/{sample}/AUGUSTUS.PASS2/{sample}.final.gff3",
		best_aed = "results/{sample}/AUGUSTUS.PASS2/{sample}.best_aed"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		echo -e "{input}" 
		echo -e "{params.aeds}"
		for aed in $(echo -e "{params.aeds}")
		do
			echo -e "$aed\t$(grep "accuracy after optimizing" results/{params.prefix}/logs/AUGUSTUS.PASS2.$aed.{params.prefix}.stdout.txt | rev | cut -d " " -f 1 | rev)"
		done > results/{params.prefix}/AUGUSTUS.PASS2/summary.tsv

		best=$(cat results/{params.prefix}/AUGUSTUS.PASS2/summary.tsv | tr . , | sort -n -k 2 | cut -f 1 | tr , . | tail -n 1)
		echo "{params.prefix}: Best training accuracy was achieved with cutoff $best"

		ln -sf $(pwd)/results/{params.prefix}/AUGUSTUS.PASS2/$best/$best.{params.prefix} {output.best_params}
		
		#mkdir {output.best_params}
		#for f in $(ls -1 results/{params.prefix}/AUGUSTUS.PASS2/$best/$best.{params.prefix}/)
		#do
		#	ln -s $(pwd)/results/{params.prefix}/AUGUSTUS.PASS2/$best/$best.{params.prefix}/$f {output.best_params}/$(echo "$f" | sed "s/^$best.//")
		#done
		
		ln -sf $(pwd)/results/{params.prefix}/AUGUSTUS.PASS2/$best/$best.augustus.gff3 {output.gff}
		echo "$best" > {output.best_aed}

		echo -e "\n$(date)\tFinished!\n"
		"""

rule snap_pass2:
	input:
		rules.merge_MAKER_PASS1.output.all_gff,
		rules.pick_augustus_training_set.output.best_aed
	params:
#		aed = config["aed"]["snap_pass2"],
		prefix = "{sample}",
		script = "bin/snap.p2.sh"
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/SNAP.PASS2.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/SNAP.PASS2.{sample}.stderr.txt"
	output:
		snap_hmm = "results/{sample}/SNAP.PASS2/{sample}.MAKER_PASS1.snap.hmm"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		export PATH="$(pwd)/bin/maker/bin:$PATH"

		#get best aed cutoff from AUGUSTUS
		aed=$(cat {input[1]})

		if [[ ! -d results/{params.prefix}/SNAP.PASS2 ]]
		then
			mkdir results/{params.prefix}/SNAP.PASS2/
		fi
		cd results/{params.prefix}/SNAP.PASS2/

		bash $basedir/{params.script} \
		{params.prefix} \
		$basedir/{input[0]} \
		$aed \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}

		echo -e "\n$(date)\tFinished!\n"
		"""

rule initiate_MAKER_PASS2:
	input:
		maker = rules.setup_maker.output.bin,
		pred_gff = rules.pick_augustus_training_set.output.gff,
		params = rules.pick_augustus_training_set.output.best_params,
		snaphmm = rules.snap_pass2.output.snap_hmm,
		MP1_ok = rules.merge_MAKER_PASS1.output,
		gmhmm = rules.genemark.output.model,
		best_aed = rules.pick_augustus_training_set.output.best_aed
	params:
		prefix = "{sample}",
		script = "bin/prepare_maker_opts_PASS2.sh",
		protein_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.protein2genome.gff",
		rm_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.repeats.gff",
		altest_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.cdna2genome.gff",
		est_gff = "results/{sample}/MAKER.PASS1/{sample}.noseq.maker.est2genome.gff"
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS2.init.{sample}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS2.init.{sample}.stderr.txt"
	output:
		ok = "results/{sample}/MAKER.PASS2/init.ok"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		export PATH="$(pwd)/bin/maker/bin:$PATH"
		aed=$(cat {input.best_aed})

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
		$basedir/{input.gmhmm} \
		$aed.{params.prefix} \
		$basedir/{input.params} \
		$basedir/{input.pred_gff} \
		$basedir/{params.rm_gff} \
		$basedir/{params.protein_gff} \
		$basedir/{params.altest_gff} \
		$basedir/{params.est_gff} \
		$(pwd)/tmp/config \
		1> $basedir/{log.stdout} 2> $basedir/{log.stderr}
		
		retVal=$(( retVal + $? ))


		if [ ! $retVal -eq 0 ]
		then
			echo "There was some error" >> $basedir/{log.stderr}
			exit $retVal
		else
			touch $basedir/{output.ok}
		fi
		
		echo -e "\n$(date)\tFinished!\n"
		"""

rule run_MAKER_PASS2:
	input:
		maker = rules.setup_maker.output.bin,
		split_ok = rules.split.output.ok,
		init_ok = rules.initiate_MAKER_PASS2.output.ok
	params:
		sub = "results/{sample}/GENOME_PARTITIONS/{unit}.fasta",
		dir = "{unit}",
		prefix = "{sample}"
	threads: config["threads"]["run_MAKER_PASS2"]
	singularity:
		"docker://chrishah/premaker-plus:18"
	log:
		stdout = "results/{sample}/logs/MAKER.PASS2.run.{sample}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/MAKER.PASS2.run.{sample}.{unit}.stderr.txt"
	output:
		sub_fasta = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.fasta",
		gff = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS2/{unit}/{sample}.{unit}.noseq.maker.gff"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"

		basedir=$(pwd)
		export PATH="$(pwd)/bin/maker/bin:$PATH"

		cd results/{params.prefix}/MAKER.PASS2/{params.dir}
		ln -s $basedir/{params.sub} {params.prefix}.{params.dir}.fasta

		AUGUSTUS_CONFIG_PATH=$basedir/results/{params.prefix}/MAKER.PASS2/tmp/config
		ln -fs $basedir/results/{params.prefix}/GENEMARK/.gm_key .

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
		echo -e "\n$(date)\tStarting ...\n"
		basedir=$(pwd)
		
		cd results/{params.prefix}/MAKER.PASS2/{params.dir}/

		bash $basedir/{params.script} {params.prefix}.{params.dir}.maker.output

		echo -e "\n$(date)\tFinished!\n"
		"""	

rule merge_MAKER_PASS2:
	input:
		lambda wildcards: expand("results/{{sample}}/MAKER.PASS2/{unit}/{{sample}}.{unit}.all.maker.gff", sample=wildcards.sample, unit=unitdict[wildcards.sample])
#		expand("results/{unit.sample}/MAKER.PASS2/{unit.unit}/{unit.sample}.{unit.unit}.all.maker.gff", unit=units.itertuples())
	params:
		prefix = "{sample}",
		script = "bin/merging.sh"
	singularity:
		"docker://chrishah/premaker-plus:18"
	output:
		all_gff = "results/{sample}/MAKER.PASS2/{sample}.all.maker.gff",
		noseq_gff = "results/{sample}/MAKER.PASS2/{sample}.noseq.maker.gff",
		proteins = "results/{sample}/MAKER.PASS2/{sample}.all.maker.proteins.fasta"
	shell:
		"""
		echo -e "\n$(date)\tStarting ...\n"
		export PATH="$(pwd)/bin/maker/bin:$PATH"
		basedir=$(pwd)
		cd results/{params.prefix}/MAKER.PASS2/

		bash $basedir/{params.script} {params.prefix}

		echo -e "\n$(date)\tFinished!\n"
		"""
