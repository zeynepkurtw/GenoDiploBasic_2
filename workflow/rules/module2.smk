#Structural Annotation
rule prodigal:
    input:
         assembly="results/Genomics/1_Assembly/2_Assemblers/{assembler}/assembly.fasta",
    output:
          gff="results/Genomics/2_Annotation/1_Structural/prodigal/{assembler}/genome.gff",
          faa="results/Genomics/2_Annotation/1_Structural/prodigal/{assembler}/genome.faa",
    conda:
         "workflow/envs/genomics.yaml"
    script:
          "workflow/scripts/genomics/annotation/prodigal.py"

#Functional Annotation
rule make_diamond_db:
    input:
        "resources/DiploProteoms/{db}.fa"
    output:
        "resources/DiploProteoms/{db}.db.dmnd"
    conda:
        "workflow/envs/genomics.yaml"
    shell:
        "diamond makedb --in {input} --db {output}"
    #script:
     #   "workflow/scripts/genomics/annotation/makedb.py"

rule diamond_blastp:
    input:
        genome="results/Genomics/2_Annotation/1_Structural/{annotation}/{assembler}/genome.faa",
        db= "resources/DiploProteoms/{db}.db.dmnd"
    output:
        "results/Genomics/2_Annotation/2_Functional/blastp/{annotation}/{assembler}/{db}/genome.blastp"
    params:
        outfmt="6 qseqid sseqid evalue qlen slen length pident stitle",
        evalue=0.00001,
        threads=32,
        max_target_seqs=1,
        max_hsps=1,
        more_sensitive="-b5 -c1"
    conda:
        "workflow/envs/genomics.yaml"
    script:
        "workflow/scripts/genomics/annotation/diamond.py"

