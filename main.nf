raw_reads = Channel.fromFilePairs(params.reads)

process trimReads {
  input:
  set val(id), file(reads) from raw_reads

  output:
  set val(id), file("*val_1.fq.gz"), file ("*val_2.fq.gz") into trimmed_reads
  file "*trimming_report.txt" into trimming_reports
  file "*fastqc*" into fastqc_reports

  publishDir "$id-output/trim_galore"
  container "picozoa.sif"
  maxForks 2

  script:
  """
  trim_galore --paired --fastqc --gzip $reads
  """
}

process assembly {
  input:
  set val(id), file(forward), file (reverse) from trimmed_reads

  output:
  set val(id), file("spades/contigs.fasta") into assembled_contigs

  publishDir "$id-output"
  cpus 20
  memory '50 GB'
  maxForks 1

  script:
  """
  mkdir tmp_dir
  /opt/SPAdes-3.13.0-Linux/bin/spades.py -o spades --threads ${task.cpus} \
                      --tmp-dir tmp_dir -1 $forward -2 $reverse \
                      --sc --careful -k 21,33,55,99 \
                      --memory 50
  """
}

assembled_contigs.into{assembled_contigs_quast; assembled_contigs_prodigal; assembled_contigs_barrnap}

process quast {
  input:
  set val(id), file(contigs) from assembled_contigs_quast

  output:
  set val(id), file("quast/*") into quast_results

  publishDir "$id-output"
  cpus 4

  script:
  """
  /opt/quast-4.5/quast.py -t ${task.cpus} -m 2000 -o quast $contigs
  """
}

process predictGenes {
  input:
  set val(id), file(contigs) from assembled_contigs_prodigal

  output:
  file("${contigs.simpleName}.gff") into prodigal_gff
  set val(id), file("${contigs.simpleName}.fna") into prodigal_fna
  set val(id), file("${contigs.simpleName}.faa") into prodigal_faa

  publishDir "$id-output/prodigal"
  container "picozoa.sif"

  script:
  """
  prodigal -f gff -o ${contigs.simpleName}.gff -d ${contigs.simpleName}.fna -a ${contigs.simpleName}.faa -p meta -i $contigs
  """
}

process predictRNA {
  input:
  set val(id), file(contigs) from assembled_contigs_barrnap

  output:
  set val(id), file("*.gff") into barrnap_predictions optional true

  publishDir "$id-output/barrnap"
  container "picozoa.sif"
  cpus 4

  script:
  """
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom bac $contigs > bac.gff
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom arc $contigs > arc.gff
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom mito $contigs > mito.gff
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom euk $contigs > euk.gff
  """
}

process diamondNR {
  input:
  set val(id), file(queries) from prodigal_fna

  output:
  set val(id), file("${queries.simpleName}.daa") into diamond_daa

  publishDir "$id-output/diamond"
  cpus 10
  maxForks 3

  script:
  """
  /opt/diamond_0.9.9/diamond blastx --threads 30 \
                -f 100 \
                --sensitive \
                -o ${queries.simpleName}.daa \
                -q $queries \
                -d /media/Data_2/diamond_dbs/nr_0.9.9.dmnd
  """
}



process DAAMeganizer {
  input:
  set val(id), file(daa) from diamond_daa

  output:
  set val(id), file(daa) into meganized_daa

  publishDir "$id-output/megan"
  cpus 10
  maxForks 3

  script:
  """
  /opt/megan/tools/daa-meganizer --in $daa \
                                 --acc2taxa /media/Data_2/megan/prot_acc2tax-Nov2018X1.abin \
                                 --lcaAlgorithm Weighted
  """
}
