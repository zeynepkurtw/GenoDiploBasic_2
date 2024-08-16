rule fastqc:
    input:
         input_dir= lambda wildcards: config["data_dir"],
    params:
          threads=32,
    output:
          out_dir=directory("results/Genomics/1_Assembly/1_Preprocessing/fastqc/"),
    conda:
         "envs/genomics.yaml",
    script:
          "scripts/Genomics/1_Assembly/1_Preprocessing/ReadQualityCheck.py"

#Assembly
rule flye:
    input:
         reads=lambda wildcards: config["data_dir"],
    params:
          #genome_size="114m",
          genome_size = "27m",
          threads=32,
    output:
           out_dir= "results/Genomics/1_Assembly/2_Assemblers/flye/assembly.fasta"
           #out_dir = directory("results/Genomics/1_Assembly/2_Assemblers/flye/"),
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/2_Assemblers/FlyeAssembler.py"

"""rule setup_nr_db:  #FIX how to actuvste this before running blastn
    input:
        #outdir = protected(directory("/data/zeynep/databases"))
        outdir = lambda wildcards: config["blast_db"]
    output:
        "results/Genomics/1_Assembly/3_Evaluation/nr_db/nr_db.dmnd"
    conda:
        "envs/genomics.yaml"
    script:
        "scripts/Genomics/1_Assembly/3_Evaluation/setup_nr_db.py"""

rule blastn:
    input:
        query="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
        #db="/data/zeynep/databases"
        db= lambda wildcards: config["blast_db"]
    output:
        "results/Genomics/1_Assembly/3_Evaluation/blastn/{assembler}/{db}/assembly.blastn"
    params:
        outfmt= "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle",
        threads=32,
        evalue=1e-10,
        db_prefix="/data/zeynep/databases/{db}"
    conda:
        "envs/genomics.yaml"
    script:
        "scripts/Genomics/1_Assembly/3_Evaluation/blastn.py"

rule meryl:
    input:
         genome="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
    output:
          merylDB=directory("results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/merlyDB"),
          repetitive_k15="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/repetitive_k15.txt",
    params:
          threads=30,
          nanopore=True
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/CalculateKmerLongReads.py"

rule winnowmap:
    input:
         genome="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
         long_read="/data/zeynep/barkhanus_data/DNA/raw/{long_read}.fastq.gz",
         merylDB="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/merlyDB",
         repetitive_k15="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/repetitive_k15.txt",
    output:
          sorted_bam="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/{long_read}.bam",
    params:
          threads=32,
          nanopore=True
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/MapLongReadsToAssembly.py"

#Evaluation
rule quast:
    input:
         assembly="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
    params:
          threads=32
    output:
          report_dir=directory("results/Genomics/1_Assembly/3_Evaluation/quast/{assembler}/")
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/AssemblyQualityCheck.py"

rule multiqc:
    input:
         input_dir="results/Genomics/1_Assembly/2_Assemblers/{assembler}/",
    params:
          threads=32
    output:
          out_dir=directory("results/Genomics/1_Assembly/3_Evaluation/multiqc/{assembler}")
    conda:
         "envs/genomics.yaml"
    shell:
         'multiqc {input.input_dir} -o {output.out_dir}'

rule plot_coverage_cont:
    input:
         #coverage on assembley
         nano="results/Genomics/1_Assembly/3_Evaluation/winnowmap/{assembler}/nanopore.bam"
    output:
          out="results/Genomics/1_Assembly/3_Evaluation/deeptools/{assembler}.png",
          outraw="results/Genomics/1_Assembly/3_Evaluation/deeptools/{assembler}/outRawCounts.txt"
    params:
          threads=32,
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/Genomics/1_Assembly/3_Evaluation/PlotCoverage.py"

