workflow mm10_workflow {
    String rna_dna_fastq_file

    #where rna_dna_fastq_file contains fastq files in the form sample_name sample1_rna_r1.fastq.gz sample2_rna_r2.fastq.gz sample1_dna_r1.fastq.gz sample2_dna_r2.fastq.gz(one row per sample)
    Array[Array[File]] rna_dna_fastqs = read_tsv(rna_dna_fastq_file)



    String star_index_dir = "/home/exacloud/gscratch/mcweeney_lab/jengs/genomes/STAR_idx/STAR_mm10"    
    String bwa_index_dir = "/home/exacloud/gscratch/mcweeney_lab/jengs/genomes/mm10"
    File ref_genome_dir = "/home/exacloud/gscratch/mcweeney_lab/jengs/genomes/mm10/GRCm38.primary_assembly.genome.fa"
    File ref_genome_gtf = "/home/exacloud/gscratch/mcweeney_lab/jengs/genomes/Strand_detection/mm10.ncbiRefSeq_sorted.gtf"
	
	scatter (fq in rna_dna_fastqs) {
		call fastp {
			input:
				cur_sample=fq[0],
				r1=fq[1],
				r2=fq[2]
		}
		
		call star{
			input:
			    	cur_sample=fq[0],
				#input_fastqs=fastp.out_fastqs,
				read1=fastp.read1,
				read2=fastp.read2,
				star_index_dir=star_index_dir
		}
                call bwa_mem {
                	input:
                        	cur_sample = fq[0],
                                r1 = fq[3],
                                r2 = fq[4],
				ref_genome_dir=ref_genome_dir
                }
                call samtools {
                	input:
                        	cur_sample = fq[0],
                        	aligned_sam = bwa_mem.out_sam
                }
		#call samtools_index {
		#	input:
		#		aligned_bam=star.out_bam
		#}

		#call reditools2_rna {
			#input:
				#cur_sample=fq[0],
				#aligned_bam=star.out_bam,
				#aligned_bam_bai=star.out_bai,
				#ref_genome_dir=ref_genome_dir
		#}
	call reditool_rna_dna {
		input:
			cur_sample=fq[0],
			input_rna=star.out_bam,
			input_rna_bai = star.out_bai,
			input_dna=samtools.out_bam,
			input_dna_bai=samtools.out_bai,
			ref_genome_dir=ref_genome_dir,
			ref_genome_gtf=ref_genome_gtf		
	}
	}
	output{
        	#Array[File] reditools2_rna_out = reditools2_rna.out_txt
		Array[File] reditools_rna_dna_out = reditool_rna_dna.out_txt
        	#Array[File] star_bam = star.out_bam
		#Array[File] star_bam_bai = star.out_bai
		#Array[File] bwa_bam = samtools.out_bam
		#Array[File] bwa_bai = samtools.out_bai
    }
}
task fastp {
	String cur_sample
        File r1
	File r2
        command {
                set -e
		fastp -i ${r1}  -I ${r2} -o \
		${cur_sample}_trimmed_r1.fastq.gz -O ${cur_sample}_trimmed_r2.fastq.gz -q 25 -u 10 -l 50 -y -x -w 4
	}
        runtime {
        	runtime_minutes: 600
        	requested_memory_per_core: "30G"
        	cpus: 9
        	maxRetries: 1
    	}   
        output {
                File read1 = "${cur_sample}" + "_trimmed_r1.fastq.gz"
		File read2 ="${cur_sample}" + "_trimmed_r2.fastq.gz"
        }       
}       

task star{
	String cur_sample
	File read1
	File read2
	String star_index_dir
	
	command {
		set -e
		STAR --runMode alignReads --runThreadN 4 --genomeDir ${star_index_dir} \
			--genomeLoad NoSharedMemory \
                        --outReadsUnmapped Fastx \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMstrandField intronMotif \
                        --outSAMattributes All  \
                        --readFilesIn ${read1} ${read2} \
                        --readFilesCommand "gunzip -c" \
                        --outFilterType BySJout \
                        --outFilterMultimapNmax 1 \
                        --alignSJoverhangMin 8 \
                        --alignSJDBoverhangMin 1 \
                        --outFilterMismatchNmax 999 \
                        --outFilterMismatchNoverLmax 0.04 \
                        --alignIntronMin 20 \
                        --alignIntronMax 1000000 \
                        --alignMatesGapMax 1000000 \
			--outFileNamePrefix "RNA_Sample_${cur_sample}_"  
			/home/exacloud/gscratch/mcweeney_lab/jengs/software/samtools/samtools-1.14/samtools index RNA_Sample_${cur_sample}_Aligned.sortedByCoord.out.bam

	}
	runtime {
        runtime_minutes: 600
        requested_memory_per_core: "30G"
        cpus: 9
        maxRetries: 1
    	}
	output {
		File out_bam = "RNA_Sample_" + "${cur_sample}" + "_Aligned.sortedByCoord.out.bam"
		File out_bai = "RNA_Sample_" + "${cur_sample}" + "_Aligned.sortedByCoord.out.bam.bai"
	}
}
task bwa_mem {
                String cur_sample
                File r1
                File r2
		String ref_genome_dir
        command {
                set -e
                /home/groups/mcweeney_lab/jengs/software/bwa-0.7.17/bwa mem ${ref_genome_dir} -Y ${r1} ${r2} > ${cur_sample}.sam
        }
        runtime {
                runtime_minutes: 600
                requested_memory_per_core: "30G"
                cpus: 9
                maxRetries: 1
        }
        output {
                File out_sam="${cur_sample}" + ".sam"
        }

}
task samtools {
                String cur_sample
                File aligned_sam
        command {
                set -e
                /home/exacloud/gscratch/mcweeney_lab/jengs/software/samtools/samtools-1.14/samtools view -S ${aligned_sam} -b -o ${cur_sample}.bam
                /home/exacloud/gscratch/mcweeney_lab/jengs/software/samtools/samtools-1.14/samtools sort ${cur_sample}.bam -o ${cur_sample}_sorted.bam -O bam -T tmp
                /home/exacloud/gscratch/mcweeney_lab/jengs/software/samtools/samtools-1.14/samtools index ${cur_sample}_sorted.bam
        }
        runtime {
                runtime_minutes: 600
                requested_memory_per_core: "30G"
                cpus: 9
                maxRetries: 1
        }
        output {
                File out_bam = "${cur_sample}" + "_sorted.bam"
                File out_bai = "${cur_sample}" + "_sorted.bam.bai"
        }

}

task reditools2_rna {
	String cur_sample
	#File aligned_bam
	#File aligned_bam_bai
	File ref_genome_dir

	command {
		set -e
		python2.7 /home/exacloud/gscratch/mcweeney_lab/jengs/software/REDItools2/src/cineca/reditools.py -f /home/exacloud/gscratch/mcweeney_lab/jengs/RNA210819JS/RNA210819JS/211015_A01058_0175_AHTWCCDSX2/RNA210819JS/rna_bam_files/RNA_Sample_${cur_sample}_Aligned.sortedByCoord.out.bam -r ${ref_genome_dir} -o RNA_Sample_${cur_sample}_redi_rna_table.txt
	}
    runtime {
        runtime_minutes: 600
        requested_memory_per_core: "10G"
        cpus: 1
        maxRetries: 1
    }
    output {
       File out_txt = "${cur_sample}" + "_redi_rna_table.txt"
    }

}
task reditool_rna_dna {
        String cur_sample
        File input_rna
        File input_rna_bai
        File input_dna
        File input_dna_bai         
	File ref_genome_dir
        File ref_genome_gtf
        command {
        	set -e
	python2.7 /home/exacloud/gscratch/mcweeney_lab/jengs/software/REDItools/main/REDItoolDnaRna.py -i /home/exacloud/gscratch/mcweeney_lab/jengs/RNA210819JS/RNA210819JS/211015_A01058_0175_AHTWCCDSX2/RNA210819JS/rna_bam_files/RNA_Sample_${cur_sample}_Aligned.sortedByCoord.out.bam -j /home/exacloud/gscratch/mcweeney_lab/jengs/RNA210819JS/RNA210819JS/211015_A01058_0175_AHTWCCDSX2/RNA210819JS/dna_bam_files//${cur_sample}_sorted.bam -o ${cur_sample} -f ${ref_genome_dir} -t10 -c1,1 -m30,255 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S -G ${ref_genome_gtf}        
        }
   	runtime {
        	runtime_minutes: 2160
        	requested_memory_per_core: "10G"
       		cpus: 10
        	maxRetries: 1
    	}

        output {
        	File out_txt = "${cur_sample}" + "_redi_rna_dna_table.txt"
        }
}
