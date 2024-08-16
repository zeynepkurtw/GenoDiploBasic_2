rule fastqc:
    input:
         input_dir= lambda wildcards: config["data_dir"],
    params:
          threads=32,
    output:
          out_dir=directory("results/Genomics/1_Assembly/1_Preprocessing/fastqc/"),
    conda:
         "envs/genomics.yaml"
    script:
          "scripts/genomics/assembly/fastqc.py"

#Assembly
rule flye:
    input:
         reads=lambda wildcards: config["data_dir"],
    params:
          #genome_size="114m",
          genome_size = "27m",
          threads=32,
    output:
           out_dir = directory("results/Genomics/1_Assembly/2_Assemblers/flye/"),
           out= "results/Genomics/1_Assembly/2_Assemblers/flye/assembly.fasta",
    conda:
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/assembly/flye.py"

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
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/assembly/meryl.py"

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
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/assembly/winnowmap.py"

#Evaluation
rule quast:
    input:
         assembly="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
    params:
          threads=32
    output:
          report_dir=directory("results/Genomics/1_Assembly/3_Evaluation/quast/{assembler}/")
    conda:
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/assembly/quast.py"

rule multiqc:
    input:
         input_dir="results/Genomics/1_Assembly/2_Assemblers/{assembler}/",
    params:
          threads=32
    output:
          out_dir=directory("results/Genomics/1_Assembly/3_Evaluation/multiqc/{assembler}")
    conda:
         "workflow/envs/genomics.yaml"
    script:
         "workflow/scripts/genomics/assembly/multiqc.py"

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
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/assembly/deeptools.py"

