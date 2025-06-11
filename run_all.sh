#!/bin/bash


# １．遺伝子名を指定してgnomadからバリアントリストを取ってくる
./scripts/01_get_gnomad_variants_by_gene.py \
  --gene_symbol RNU4-2 \
  --reference_genome GRCh38 \
  --dataset gnomad_r4 \
  --out_vcf results/gnomad_rnu4-2.vcf


# ２．そのリストを与えつつ、ゲノム領域を指定してさらに絞って、gnomadバリアントの分布を絵にする
./scripts/02_plot_gnomad_variants.R \
  --vcf results/gnomad_rnu4-2.vcf \
  --ref_ver GRCh38 \
  --chr chr12 \
  --start 120291763 \
  --end 120291903 \
  --min_ac 1 \
  --out results/gnomad_rnu4-2.png


# ３．ClinVarのvcfを与えて、ゲノム領域を指定して、バリアントの分布を絵にする
# ついでに、ドメインの領域のファイルを与えて、絵にする
./scripts/03_plot_clinvar_variants.R \
  --clinvar_vcf data/clinvar.vcf.gz \
  --reference_genome GRCh38 \
  --chr 12 \
  --start 120291763 \
  --end 120291903 \
  --domain_file data/rnu4-2_domain.txt \
  --out_png results/clinvar_rnu4-2.png
