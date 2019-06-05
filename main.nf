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
  memory '30 GB'
  maxForks 3

  script:
  """
  mkdir tmp_dir
  /opt/SPAdes-3.13.0-Linux/bin/spades.py -o spades --threads ${task.cpus} \
                      --tmp-dir tmp_dir -1 $forward -2 $reverse \
                      --sc --careful -k 21,33,55,99 \
                      --memory 50
  """
}

assembled_contigs.into{assembled_contigs_quast; assembled_contigs_prodigal; assembled_contigs_barrnap; assembled_contigs_stats}

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
  set val(id), file("bac.gff"), file(contigs) into barrnap_predictions_bac optional true
  set val(id), file("euk.gff"), file(contigs) into barrnap_predictions_euk optional true

  publishDir "$id-output/barrnap"
  container "picozoa.sif"
  cpus 4

  script:
  """
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom bac $contigs > bac.gff
  barrnap --threads ${task.cpus} --reject 0.2 --kingdom euk $contigs > euk.gff
  for gff in *.gff; do if grep -q -v '^#' \$gff; then echo \$gff; else rm \$gff; fi ; done
  """
}

barrnap_predictions_bac.concat( barrnap_predictions_euk ).into{barrnap_predictions_SSU; barrnap_predictions_LSU}

process classifySSU {
  input:
  set val(id), file(gff), file(contigs) from barrnap_predictions_SSU

  output:
  set val(id), file("${gff.simpleName}.ssu.lca") into rna_lca optional true
  set val(id), file("${gff.simpleName}.ssu.blastn") into rna_blastn optional true
  set val(id), file("${gff.simpleName}.ssu.fna") into rna_fasta optional true

  publishDir "$id-output/barrnap"
  cpus 4


  script:
  """
  python3 /media/Data_1/Max/misc-scripts/parse_barrnap.py -t ${task.cpus} -l $gff -r $contigs -f 'SSU' -o ${gff.simpleName}.ssu.fna
  if [ -e ${gff.simpleName}.ssu.fna ];
  then
      /opt/ncbi-blast-2.7.1+/bin/blastn -db /media/Data_2/SILVA/v132/SILVA_132_SSURef_Nr99_tax_silva_trunc.db \
                    -query ${gff.simpleName}.ssu.fna \
                    -num_threads ${task.cpus} \
                    -out ${gff.simpleName}.ssu.blastn

      /opt/megan/tools/blast2lca -i ${gff.simpleName}.ssu.blastn \
                    -f BlastText \
                    -m BlastN \
                    -o ${gff.simpleName}.ssu.lca \
                    -s2t /media/Data_2/megan/SSURef_Nr99_132_tax_silva_to_NCBI_synonyms.map.gz
  fi
  """
}


process classifyLSU {
  input:
  set val(id), file(gff), file(contigs) from barrnap_predictions_LSU

  output:
  set val(id), file("${gff.simpleName}.lsu.lca") into rna_lca_lsu optional true
  set val(id), file("${gff.simpleName}.lsu.blastn") into rna_blastn_lsu optional true
  set val(id), file("${gff.simpleName}.lsu.fna") into rna_fasta_lsu optional true

  publishDir "$id-output/barrnap"
  cpus 4


  script:
  """
  python3 /media/Data_1/Max/misc-scripts/parse_barrnap.py -t ${task.cpus} -l $gff -r $contigs -f 'LSU' -o ${gff.simpleName}.lsu.fna
  if [ -e ${gff.simpleName}.lsu.fna ];
  then
      /opt/ncbi-blast-2.7.1+/bin/blastn -db /media/Data_2/SILVA/v132/SILVA_132_LSURef_tax_silva_trunc.db \
                    -query ${gff.simpleName}.lsu.fna \
                    -num_threads ${task.cpus} \
                    -out ${gff.simpleName}.lsu.blastn

      /opt/megan/tools/blast2lca -i ${gff.simpleName}.lsu.blastn \
                    -f BlastText \
                    -m BlastN \
                    -o ${gff.simpleName}.lsu.lca \
                    -s2t /media/Data_2/megan/LSURef_132_tax_silva_to_NCBI_synonyms.map.gz
  fi
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
  /opt/diamond_0.9.9/diamond blastx --threads ${task.cpus} \
                -f 100 \
                --sensitive \
                -o ${queries.simpleName}.daa \
                -q $queries \
                -d /media/Data_2/diamond_dbs/nr_0.9.9.dmnd
  """
}

diamond_daa.into{diamond_daa_megan; diamond_daa_blast}

process diamondToBlast {
  input:
  set val(id), file(daa) from diamond_daa_blast

  output:
  set val(id), file("${daa.simpleName}.tab") into diamond_tab

  publishDir "$id-output/diamond"
  cpus 10
  maxForks 3

  script:
  """
  /opt/diamond_0.9.9/diamond view --daa ${daa} \
                --threads ${task.cpus} \
                -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                -o ${daa.simpleName}.tab \
  """
}


process DAAMeganizer {
  input:
  set val(id), file(daa) from diamond_daa_megan

  output:
  set val(id), file(daa) into meganized_daa

  stageInMode 'copy'
  publishDir "$id-output/megan-adjusted"
  cpus 10
  maxForks 3

  script:
  """
  /opt/megan/tools/daa-meganizer --in $daa \
                                 --acc2taxa /media/Data_2/megan/prot_acc2tax-Nov2018X1.abin \
                                 --lcaAlgorithm Weighted \
                                 --maxExpected 0.001 \
                                 --minPercentIdentity 80.0 \
                                 --topPercent 40
  """
}


process DaaToInfo {
  input:
  set val(id), file(daa) from meganized_daa

  output:
  set val(id), file("r2c.tab") into r2c
  set val(id), file("c2c.tab") into c2c

  publishDir "$id-output/megan-adjusted"

  script:
  """
  /opt/megan/tools/daa2info -r2c Taxonomy \
                            -i $daa >> r2c.tab

  /opt/megan/tools/daa2info -c2c Taxonomy \
                            -i $daa >> c2c.tab
  """
}

tax_names = r2c.combine(assembled_contigs_stats)


process getTaxNames {
  input:
  set val(id), file(r2c), val(id2), file(contigs) from tax_names

  output:
  set val(id), file("contig2tax.tab") into contig2tax

  publishDir "$id-output/megan-adjusted"

  when:
  "${id}" == "${id2}"

  script:
  """
  #! /usr/bin/python3

  import ete3
  from ete3 import ncbi_taxonomy
  ncbi = ncbi_taxonomy.NCBITaxa()
  from collections import defaultdict
  from Bio import SeqIO

  name2tax = defaultdict(list)
  name2group = defaultdict(set)
  name2length = defaultdict(str)
  groups = ["Eukaryota", "Bacteria", "Archaea", "Viruses"]
  for line in open("${r2c}"):
      line = line.strip().split()
      name = "_".join(line[0].split('_')[:-1])
      taxid = int(line[1])
      try:
          tax = ncbi.get_taxid_translator([taxid])[taxid]
          lineage = ncbi.get_taxid_translator(ncbi.get_lineage(taxid)).values()
          for g in groups:
              if g in lineage:
                  name2group[name].add(g)
          name2tax[name].append(tax)
          name2length[name] = line[0].split('_')[3]
      except:
          print("Could not find taxonomy")

  for rec in SeqIO.parse("$contigs", 'fasta'):
      if not rec.id in name2tax.keys():
          name2length[rec.id] = rec.id.split('_')[3]
          name2tax[rec.id] = ["no hit"]
          name2group[rec.id] = set(["no hit"])

  with open('contig2tax.tab', 'w') as out:
      for k, v in name2tax.items():
          print("\\t".join([k, name2length[k], ";".join(name2group[k]), ";".join(v)]), file=out)
  """
}
