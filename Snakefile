# Setting up the variables that will be used
snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar"
snpeff_genome = 'hg38'
snpeff_db_folder = '/staging/leuven/stg_00079/teaching/snpeff_db'
genome_db = "/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa"

# Fetch all files and strip the .fq.gz extension from the file names
sample_names_with_ext, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")
filtered_samples = [name.rsplit('.fq.gz', 1)[0] for name in sample_names_with_ext if re.match(r"HG0\d{3}9", name)]

# Print for check the files that are going to be analyzed
print(filtered_samples)

# ALL - the first rule is the rule snakemake automatically executes
rule all:
    input:
        raw_copied=expand("000.base/{file}.fq.gz", file=filtered_samples),
        fastqc_zip=expand("010.fastqc/{file}_fastqc.zip", file=filtered_samples),
        rep1=expand("010.fastqc/{file}_fastqc/Images/per_base_quality.png", file=filtered_samples),
        vcf="030.samtools/snps.vcf",
        cleaned_vcf="040.cleaned/snps.cleaned.vcf",
        snpeff ="050.snpeff/snps.annotated.vcf",
        genes = "genes.vcf",
        report="reports/workflow_report.log",
        rep4="040.cleaned/cleaned_quality_scores_histogram.png",
        rep5="reports/gene_snp_heatmap.png",
        

# Create the report with the print statements after the execution of a rule
rule init_report:
    output:
        report=report("reports/workflow_report.log", category="Workflow_Report"),
    shell:
        "touch {output.report}"

# Copying the sample files to a new directory
rule copy_files:
    input:
        raw=expand("/staging/leuven/stg_00079/teaching/1000genomes/{file}.fq.gz", file=filtered_samples),
        report_init="reports/workflow_report.log",
    output:
        raw_copied="000.base/{file}.fq.gz",
    shell:
        """
        echo "Copying file: {input.raw}"
        mkdir -p 000.base
        cp {input.raw} 000.base

        # Test if the copied file exists or if they are empty and if yes report it.
        if [ ! -s {output.raw_copied} ]; then
            echo "ERROR: Copy of {wildcards.file} does not exist or is zero bytes." >> {input.report_init}
        else
            echo "OK: Copy of {wildcards.file} exists and it's not zero bytes." >> {input.report_init}
        fi
        """

# Check the quality of the samples and find any failures. Also include some diagnostic images in the report
rule fastqc:
    input:
        fq=expand("000.base/{file}.fq.gz", file=filtered_samples),
        report_init="reports/workflow_report.log",
    output:
        fastqc_zip="010.fastqc/{file}_fastqc.zip",
        html="010.fastqc/{file}_fastqc.html",
        summarydata="010.fastqc/{file}_fastqc/fastqc_data.txt",
        rep1=report("010.fastqc/{file}_fastqc/Images/per_base_quality.png", category="Fastqc",
                    subcategory="Per base quality", labels={"sample": "{file}"}),
        rep2=report("010.fastqc/{file}_fastqc/Images/per_base_sequence_content.png", category="Fastqc",
                     subcategory="Per base sequence content", labels={"sample": "{file}"}),
        rep3=report("010.fastqc/{file}_fastqc/summary.txt", category="Fastqc",
                    subcategory="Summary text", labels={"sample": "{file}"}),
    shell:
        """
        echo "Input Fastq: {input.fq} "
        fastqc -o 010.fastqc {input.fq} --extract

        # Check for FAIL statuses in the FastQC output and report them.
        if grep -q 'FAIL' {output.rep3}; then
            echo "Quality check FAILED for {wildcards.file}" >> {input.report_init}
        else
            echo "Quality check passed for {wildcards.file}" >> {input.report_init}
        fi
        """

# Mapping the reads with the database
rule bwa:
    input:
        fq=lambda wildcards: "000.base/{}.fq.gz".format(wildcards.file),
        report_init="reports/workflow_report.log",
    output:
        bam="020.bwa/{file}.bam",
        bai="020.bwa/{file}.bam.bai"
    params:
        db=genome_db
    shell:
        """
        bwa mem {params.db} {input.fq} | samtools sort -o {output.bam}
        samtools index {output.bam}

        # Test if the bam file exists or if it's empty and if yes report it 
        if [ ! -s {output.bam} ]; then
            echo "ERROR: BAM file for {wildcards.file} does not exist or is zero bytes." >> {input.report_init}
        else
            echo "OK: BAM file for {wildcards.file} exists and it's not zero bytes." >> {input.report_init}
        fi
        """

# Calling genetic variants from the aligned reads
rule variant_calling:
    input:
        db=genome_db,
        bams=expand("020.bwa/{file}.bam", file=filtered_samples),
        report_init="reports/workflow_report.log",
    output:
        vcf="030.samtools/snps.vcf",
    shell:
        """
        bcftools mpileup -Ou -f {input.db} {input.bams} \
             | bcftools call -mv -Ov -o {output.vcf}

         # Check if the VCF file exists and is greater than zero bytes, if not, report it
        if [ ! -s {output.vcf} ]; then
            echo "ERROR: VCF file does not exist or is zero bytes." >> {input.report_init}
        else
            echo "OK: VCF file exists and it is not zero bytes." >> {input.report_init}
        fi

        # Check that all variants are located in chromosome 21, if not, report it
        awk '!/^#/ {{if ($1 != "chr21") {{print "ERROR: Found SNP not in chr21 at line: " NR;}}}}' {output.vcf}
        if [ $? -ne 0 ]; then
            echo "ERROR: Some variants are not in chromosome 21." >> {input.report_init}
        else
            echo "OK: All variants are in chromosome 21." >> {input.report_init}
        fi
        """
# Keeping high-confidence variants 
rule variant_cleanup:
    input:
        db=genome_db,
        vcf="030.samtools/snps.vcf",
        report_init="reports/workflow_report.log",
    output:
        vcf="040.cleaned/snps.cleaned.vcf"
    shell:
        r"""
        # Decompose, normalize, and filter variants
        cat {input.vcf} \
            | vt decompose - \
            | vt normalize -n -r {input.db} - \
            | vt uniq - \
            | vt view -f 'QUAL>20' -h - \
            > {output.vcf}

        # Verify that all variants have QUAL>=20 and print violations if any and report them
        awk '!/^#/ && $6 < 20 {{ print "ERROR: Found variant with QUAL<20: " $0;}}' {output.vcf}
        if [ $? -ne 0 ]; then
            echo "ERROR: One or more variants are below quality threshold detected." >> {input.report_init}
        else
            echo "OK: All variants are above quality threshold." >> {input.report_init}
        fi
        """
# Annotating the cleaned variants
rule snpeff:
    input:
        vcf = "040.cleaned/snps.cleaned.vcf",
        report_init="reports/workflow_report.log",
    params:
        snpeff_db_folder = snpeff_db_folder,
        snpeff_jar = snpeff_jar,
    log:
        err="050.snpeff/snakemake.err",
    output:
        vcf = "050.snpeff/snps.annotated.vcf",
        html = "050.snpeff/snpEff_summary.html",
        genetxt = "050.snpeff/snpEff_genes.txt",
    shell:
        """

        mkdir -p 050.snpeff

        java -Xmx4096m -jar \
            {params.snpeff_jar} eff hg19 \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.vcf}

        # move output files to the snpeff output folder
        mv snpEff_genes.txt snpEff_summary.html 050.snpeff

        # Check for the presence of 'ANN' attribute in the VCF and report it
        if ! grep -q '##INFO=<ID=ANN' {output.vcf}; then
            echo "ERROR: VCF file does not contain expected SnpEff annotations." >> {input.report_init}
        else
            echo "OK: VCF file contains expected SnpEff annotations." >> {input.report_init}
        fi
        """
# Extracting the variants that are only inside the three genes of interest (intergenic regions are excluded).
rule extract_genes_of_interest:
    input:
        vcf="050.snpeff/snps.annotated.vcf",
        report_init="reports/workflow_report.log",
    output:
        vcf="genes.vcf",
    shell:
        """
        #grep -E '^#|APP|SOD1|DYRK1A' {input.vcf} > {output.vcf}
        #The "gene field" with the genes of interest is located at the "info column"
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
                        if (fields[4] == "APP" || fields[4] == "SOD1" || fields[4] == "DYRK1A") {{ 
                            print; 
                            break; 
                        }}
                    }}
                }} 
            }}
        }}' {input.vcf} > {output.vcf}

        # Check that only expected genes are in the output VCF and report
        if ! awk -F'\t' 'BEGIN {{ found=0 }} 
        !/^#/ && (index($0, "APP") == 0 && index($0, "SOD1") == 0 && index($0, "DYRK1A") == 0) {{ 
            found=1; 
            print "Unexpected line in output:", $0;  
        }}
        END {{ if (found) print "Validation failed: Non-target genes found in the output." }}' {output.vcf}; then
            echo "Check complete with issues: Non-target genes found." >> {input.report_init}
        else
            echo "Validation successful: All genes in the output are expected." >> {input.report_init}
        fi
        """

# Rule as a extra check for the cleaned variants. The histogram that is generated shows that all variants has quality>20 and it's in the report
rule plot_cleaned_quality_scores:
    input:
        cleaned_vcf="040.cleaned/snps.cleaned.vcf",
    output:
        rep4=report("040.cleaned/cleaned_quality_scores_histogram.png", category = "cleaned snps"),
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        def extract_quality_scores(vcf_path):
            scores = []
            with open(vcf_path, 'r') as file:
                for line in file:
                    if not line.startswith('#'):
                        parts = line.split('\t')
                        qual = float(parts[5]) if parts[5] != '.' else 0
                        scores.append(qual)
            return scores

        # Read quality scores from the cleaned VCF file
        cleaned_scores = extract_quality_scores(input['cleaned_vcf'])

        # Create a DataFrame
        df = pd.DataFrame({
            'Quality Scores': cleaned_scores
        })

        # Create the histogram plot
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Quality Scores'], kde=False, color='blue', binwidth=1)
        plt.axvline(x=20, color='red', linestyle='--', label='Quality Threshold: 20')
        plt.title('Histogram of SNP Quality Scores After Cleanup')
        plt.xlabel('Quality Scores')
        plt.ylabel('Number of SNPs')
        plt.legend()
        plt.grid(True)

        # Save the plot
        plt.savefig(output['rep4'])
        plt.close()

# Generating the bonus heatmap for each gene of interest per sample.
rule plot_gene_snp_heatmap:
    input:
        vcf="genes.vcf",
    output:
        rep5=report("reports/gene_snp_heatmap.png", category= "Heatmap"),
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Function to extract SNP counts for the specified genes per individual
        def extract_snp_counts(vcf_path, genes_of_interest):
            data = []
            with open(vcf_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    info_field = parts[7]
                    sample_data = parts[9:]
                    
                    # Extract all gene names from the ANN field in the INFO column
                    annotations = info_field.split('ANN=')[1].split(',')
                    associated_genes = set(ann.split('|')[3] for ann in annotations if ann.split('|')[3] in genes_of_interest)
                    
                    # Count SNPs for the genes of interest per sample
                    for i, sample in enumerate(sample_data):
                        genotype = sample.split(':')[0]
                        if genotype not in ['0/0', './.']:  # Non-reference and called
                            for gene in associated_genes:
                                data.append({'Individual': f'Individual {i+1}', 'Gene': gene})
            
            return pd.DataFrame(data)

        # Genes of interest
        genes_of_interest = {'APP', 'SOD1', 'DYRK1A'}

        # Process the VCF to extract SNP counts
        df = extract_snp_counts(input.vcf, genes_of_interest)

        # Pivot the data for heatmap plotting
        pivot_df = df.pivot_table(index='Gene', columns='Individual', aggfunc=len, fill_value=0)

        # Plot the heatmap
        plt.figure(figsize=(10, 4))
        sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='d')
        plt.title('SNP Association Per Gene and Individual')
        plt.ylabel('Gene')
        plt.xlabel('Individual')
        plt.tight_layout()

        # Save the plot
        plt.savefig(output.rep5)
        plt.close()
