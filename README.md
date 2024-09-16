Details

This is a snakemake workflow taking as initial input sequenced genomes by illumina technologies and giving back:
1) a .vcf file with all the SNPs related to the given genes (here related with AD)
2) a report in HTML format where it can be seen if there was any failures during the run of the workflow
3) a compressed file which contains the report in txt format and the a heatmap which associates the number of SNPs that were found to the patients.

Description

The rules that are executed during in the analysis are given with the exact order:
1) ```rule init_report``` it initiates the report that will be given at the end
2) ```rule copy_files``` it copies locally the selected raw sequenced data from the patients
3) ```rule fastqc``` uses fastqc for quality control of each sample and provides helpful images at the report
4) ```rule bwa maps``` the raw reads (using bwa samtools) with the selected database
5) ```rule variant_calling``` using bcftools mpileup and bcftools call calls the variants from the aligned reads and produces one snps.vcf file with all the unfiltered snps.
6) ```rule variant_cleanup``` using vt this rule filters and cleans all the snps from the previous step. Specifically it decomposes, normalizes and keeps the unique variants with quality score > 20.
7) ```rule snpeff``` using the imported snpeff and snpeff database the cleaned variants are annotated producing a snps.annotated.vcf, snpEff_summary.html, snpEff_genes.txt files.
8) ```rule extract_genes_of_interest``` extracts the SNPs that are related to the given genes excluding the intergenic regions, giving back the genes.vcf file.
9) ```rule plot_cleaned_quality_scores``` works as a bonus step, checking that all snps have quality above 20 and produces in the report a histogram
10) ```rule plot_gene_snp_heatmap``` takes as input the genes.vcf file and extracts the SNP counts for the specified genes per individual and plots the final heatmap.

Usage


Here are all the modifications that are required to personalize the workflow.

First of all setting up/adjusting the general variables are required without changing the name of the variable:
```
snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar" # the path to snpEff
snpeff_genome = 'hg38' # the genome that the annotoation will be based on
snpeff_db_folder = '/staging/leuven/stg_00079/teaching/snpeff_db' # the location of the database that will be used by snpEff
genome_db = "/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa" # the path to the genome that will be used for mapping the reads and variant calling
```
Then fetching all the sample files that we be analyzed are required:
```
sample_names_with_ext, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")  # here we get all the raw samples that are located to a specific directory
filtered_samples = [name.rsplit('.fq.gz', 1)[0] for name in sample_names_with_ext if re.match(r"HG0\d{3}9", name)] # from all samples that we got before here we chose to analyze only a specific proportion taking out the extention of each sample for flexibility reasons.
```
Here change the path of the selected samples:
```
rule copy_files:
    input:
        raw=expand("/staging/leuven/stg_00079/teaching/1000genomes/{file}.fq.gz", file=filtered_samples), #modify this line
        report_init="reports/workflow_report.log",

Comment this check or change the chromosome number if you know where your genes of interest are located:
Check that all variants are located in chromosome 21, if not, report it
        awk '!/^#/ {{if ($1 != "chr21") {{print "ERROR: Found SNP not in chr21 at line: " NR;}}}}' {output.vcf}
        if [ $? -ne 0 ]; then
            echo "ERROR: Some variants are not in chromosome 21." >> {input.report_init}
        else
            echo "OK: All variants are in chromosome 21." >> {input.report_init}
        fi
```
Here change the names of the genes of your interest:
```
awk 'BEGIN {{ FS=OFS="\t" }} 
        /^#/ {{ print; next }} 
        {{ 
            n = split($8, info, ";"); 
            for (i=1; i<=n; i++) {{ 
                if (info[i] ~ /^ANN=/) {{ 
                    split(info[i], ann, "="); 
                    m = split(ann[2], ann_details, ","); 
                    for (j=1; j<=m; j++) {{ 
                        split(ann_details[j], fields, "|"); 
                        if (fields[4] == "APP" || fields[4] == "SOD1" || fields[4] == "DYRK1A") {{  #modify here
                            print; 
                            break; 
                        }}
                    }}
                }} 
            }}
        }}' {input.vcf} > {output.vcf}

        # Check that only expected genes are in the output VCF and report
        if ! awk -F'\t' 'BEGIN {{ found=0 }} 
        !/^#/ && (index($0, "APP") == 0 && index($0, "SOD1") == 0 && index($0, "DYRK1A") == 0) {{ #modify here
            found=1; 
            print "Unexpected line in output:", $0;  
        }}
```
In the last rule change the genes names:
```
genes_of_interest = {'APP', 'SOD1', 'DYRK1A'}
```
To run the analysis workflow it is recommended to do it with an interractive session or even better as a batch job. The context of the slurm file as an example is in the slurm.txt file in the repository.
Changing the path to the environment with the depedencies, the parameters and the paths is required.
To submit the job use the command ```sbatch my_slurm_job.sh```

Depedencies

bwa 0.7.17-r1188

FastQC v0.12.1

bcftools 1.19

vt

SnpEff 5.2

Python 3.10.13

snakemake 7.32.4

