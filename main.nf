params.reads = ""

raw_reads = Channel.fromFilePairs(params.reads)

process trimReads {
  input:
  set val(id), file(reads) from raw_reads

  output:
  set val(id), file("*val_1.fq.gz"), file ("*val_2.fq.gz") into trimmed_reads

  publishDir "$id-output"

  script:
  """
  trim_galore --paired --fastqc --gzip $reads
  """
}

process assembly {
  input:
  set val(id), file(forward), file (reverse) from trimmed_reads

  output:

  publishDir "$id-output"
  cpus 20

  script:
  """
  mkdir tmp_dir
  spades.py -o spades --threads ${task.cpus} \
                      --tmp-dir tmp_dir -1 $forward -2 $reverse \
                      --sc --careful -k 33,55,99 \
                      --memory 50
  """
}
//
// process quast {
//
//
//
// }
//
// process predictGenes {
//
//
//
// }
//
// process predictRNA {
//
//
//
// }
