# snakemake -r -s Snakefile.py
# snakemake -s Snakefile.py --dag | dot -Tpdf > dag.pdf

import os
import yaml
import time


############################ Define some functions. ###############################
def createDir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    
    return "%d:%02d:%02d" % (hours, minutes, seconds)


def get_sample_name(path, split):
    files = os.listdir(path)
    samples = sorted(list(set([f.split(split)[0] for f in files])))
    return samples


def modify_file_extension(file_path, new_extension):
    file_name, file_extension = os.path.splitext(file_path)
    new_file_name = file_name + new_extension
    return new_file_name



########################### Read config files. ###############################
if len(config) == 0:
    if os.path.isfile("./configs/config_main.yaml"):
        configfile: "./configs/config_main.yaml"
    else:
        sys.exit("".join(["Make sure there is a config_main.yaml file in ", os.getcwd()+"/configs/", 
            " or specify one with the --configfile commandline parameter."]))
# print(config)


WORKDIR = os.getcwd()
data_path = os.path.join("data", config["fastq_dir"])
samples = get_sample_name(data_path, split = "_")

reference_genome_filename = os.path.splitext(os.path.basename(config["reference_genome_download_url"]))[0]  # remove .gz suffix
reference_annotation_filename = os.path.splitext(os.path.basename(config["reference_annotation_download_url"]))[0]  # remove .gz suffix

splice_sites_filename = modify_file_extension(reference_annotation_filename, ".ss")
exon_filename = modify_file_extension(reference_annotation_filename, ".exon")




############################# Create directories.  ############################
# start_time = time.time()
reference_genome_path = os.path.join("data", config["reference_dir"], config["reference_genome_dir"])
createDir(reference_genome_path)

reference_transcriptome_path = os.path.join("data", config["reference_dir"], config["reference_transcriptome_dir"])
createDir(reference_transcriptome_path)

reference_annotation_path = os.path.join("data", config["reference_dir"], config["reference_annotation_dir"])
createDir(reference_annotation_path)

createDir(config["log_dir"])

fastqc_output_path = os.path.join(config["output_dir"], config["output_fastqc_dir"])
createDir(fastqc_output_path)

index_output_path = os.path.join(config["output_dir"], config["output_index_dir"])
createDir(index_output_path)

alignment_output_path = os.path.join(config["output_dir"], config["output_alignment_dir"])
createDir(alignment_output_path)

sam2bam_output_path = os.path.join(config["output_dir"], config["output_sam2bam_dir"])
createDir(sam2bam_output_path)

featureCounts_output_path = os.path.join(config["output_dir"], config["output_featureCounts_dir"])
createDir(featureCounts_output_path)

# end_time = time.time()
# print("Time of creating directories: " + spend_time(start_time, end_time) + "\n\n")


############################# Download reference files. ##################################

os.system("nohup wget -c -P " + reference_genome_path + " " + config["reference_genome_download_url"] + " > logs/log_download_ref_genome.txt &")
os.system("nohup wget -c -P " + reference_annotation_path + " " + config["reference_annotation_download_url"] + " > logs/log_download_ref_annotation.txt &")
os.system("cd " + reference_genome_path + " && " + "gunzip " + os.path.basename(config["reference_genome_download_url"]))
os.system("cd " + reference_annotation_path + " && " + "gunzip " + os.path.basename(config["reference_annotation_download_url"]))
os.system("cd " + WORKDIR)



######################################### Snakemake rules #############################################
rule all:
    input:
        count = expand(os.path.join(featureCounts_output_path, "{sample_id}_count.tsv"), sample_id=samples), 
        count_summary = expand(os.path.join(featureCounts_output_path, "{sample_id}_count.tsv.summary"), sample_id=samples)


rule fastqc:
    input:
        fastq1 = expand(os.path.join(data_path, "{sample_id}_R1.fastq.gz"), sample_id=samples), 
        fastq2 = expand(os.path.join(data_path, "{sample_id}_R2.fastq.gz"), sample_id=samples)
    output:
        fastqc_1_result_html = expand(os.path.join(fastqc_output_path, "{sample_id}_R1_fastqc.html"), sample_id=samples), 
        fastqc_1_result_zip = expand(os.path.join(fastqc_output_path, "{sample_id}_R1_fastqc.zip"), sample_id=samples), 
        fastqc_2_result_html = expand(os.path.join(fastqc_output_path, "{sample_id}_R2_fastqc.html"), sample_id=samples), 
        fastqc_2_result_zip = expand(os.path.join(fastqc_output_path, "{sample_id}_R2_fastqc.zip"), sample_id=samples)
    params:
        path = fastqc_output_path
    log: 
        expand("logs/fastqc_{sample_id}.log", sample_id=samples)
    shell:
        "fastqc -t {config[num_threads]} -o {params.path} {input.fastq1} {input.fastq2}"


if config["TRIMMED"]:
    fastqc_trimmed_output_path = os.path.join(config["output_dir"], "fastqc_trimmed")
    createDir(fastqc_trimmed_output_path)

    rule trimming:
        input:
            fastq1 = os.path.join(data_path, "{sample_id}_R1.fastq.gz"), 
            fastq2 = os.path.join(data_path, "{sample_id}_R2.fastq.gz")
        output:
            fastq1_trimmed = os.path.join(fastqc_trimmed_output_path, "{sample_id}_R1_val_1.fq.gz"), 
            fastq2_trimmed = os.path.join(fastqc_trimmed_output_path, "{sample_id}_R2_val_2.fq.gz")
        params:
            outputpath = fastqc_trimmed_output_path
        log: 
            "logs/fastqc_trimmed_{sample_id}.log"
        shell:
            "trim_galore --fastqc -j {config[num_threads]} --paired --basename {wildcards.sample_id} -o {params.outputpath} {input.fastq1} {input.fastq2}"





rule indexGenome:
    input:
        reference = os.path.join(reference_genome_path, reference_genome_filename), 
        annotation_gtf = os.path.join(reference_annotation_path, reference_annotation_filename)
    output:
        splice_sites = os.path.join(index_output_path, splice_sites_filename), 
        exon = os.path.join(index_output_path, exon_filename), 
        index = os.path.join(index_output_path, modify_file_extension(reference_annotation_filename, "_tran"))
    shell:
        "hisat2_extract_splice_sites.py {input.annotation_gtf} > {output.splice_sites}"
        "&& hisat2_extract_exons.py {input.annotation_gtf} > {output.exon}"
        "&& hisat2-build --ss {output.splice_sites} --exon {output.exon} {input.reference} {output.index}"




rule alignment:
    input:
        exon = os.path.join(index_output_path, exon_filename), 
        fastq1_trimmed = os.path.join(fastqc_trimmed_output_path, "{sample_id}_R1_val_1.fq.gz"), 
        fastq2_trimmed = os.path.join(fastqc_trimmed_output_path, "{sample_id}_R2_val_2.fq.gz")
    output:
        sam = os.path.join(alignment_output_path, "{sample_id}.sam"),
        bam = os.path.join(sam2bam_output_path, "{sample_id}.bam"), 
        sorted_bam = os.path.join(sam2bam_output_path, "{sample_id}.sort.bam")
    params:
        index = os.path.join(index_output_path, modify_file_extension(reference_annotation_filename, "_tran"))
    run:
        shell("hisat2 -p {config[num_threads]} --dta -x {params.index} -1 {input.fastq1_trimmed} -2 {input.fastq2_trimmed} -S {output.sam}")
        shell("samtools view -@ {config[num_threads]} -b -S {output.sam} > {output.bam}")
        shell("samtools sort -@ {config[num_threads]} {output.bam} -o {output.sorted_bam}")


rule featureCounts:
    input:
        sorted_bam = os.path.join(sam2bam_output_path, "{sample_id}.sort.bam"), 
        annotation = os.path.join(reference_annotation_path, reference_annotation_filename)
    output:
        count = os.path.join(featureCounts_output_path, "{sample_id}_count.tsv"), 
        count_summary = os.path.join(featureCounts_output_path, "{sample_id}_count.tsv.summary")
    run:
        shell("featureCounts -p -T {config[num_threads]} -t exon -g gene_id -a {input.annotation} -o {output.count} {input.sorted_bam}")


onstart:
    print("RNA-seq workflow start.")


onsuccess:
    print("RNA-seq workflow finished, no error.")


onerror:
    print("Errors occurred!")