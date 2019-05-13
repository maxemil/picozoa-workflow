raw_reads = Channel.fromFilePairs(params.reads)

process trimReads {
  input:
  set val(id), file(reads) from raw_reads

  output:
  set val(id), file("*val_1.fq.gz"), file ("*val_2.fq.gz") into trimmed_reads

  publishDir "$id-output"
  container "picozoa.sif"

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

  script:
  """
  mkdir tmp_dir
  /opt/SPAdes-3.13.0-Linux/bin/spades.py -o spades --threads ${task.cpus} \
                      --tmp-dir tmp_dir -1 $forward -2 $reverse \
                      --sc --careful -k 33,55,99 \
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
//
// process predictGenes {
//
//
//
// }
//
process predictRNA {
  input:
  set val(id), file(contigs) from assembled_contigs_barrnap

  output:
  set val(id), file("*.gff") into barrnap_predictions optional true

  publishDir "$id-output/barrnap"
  cpus 4

  script:
  """
  /opt/barrnap-0.8/bin/barrnap --threads ${task.cpus} --reject 0.2 --kingdom bac $contigs > bac.gff
  /opt/barrnap-0.8/bin/barrnap --threads ${task.cpus} --reject 0.2 --kingdom arc $contigs > arc.gff
  /opt/barrnap-0.8/bin/barrnap --threads ${task.cpus} --reject 0.2 --kingdom mito $contigs > mito.gff
  /opt/barrnap-0.8/bin/barrnap --threads ${task.cpus} --reject 0.2 --kingdom euk $contigs > euk.gff
  """
}
