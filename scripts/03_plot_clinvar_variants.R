#!/usr/bin/env Rscript

library(VariantAnnotation)
library(trackViewer)
library(dplyr)
library(argparser)

p <- arg_parser("")
p <- add_argument(p, "--clinvar_vcf",      default='clinvar.vcf.gz',    type="character", help="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gzとtbi" )
p <- add_argument(p, "--reference_genome", default='hg38',              type="character", help="VariantAnnotationパッケージの挙動に影響するらしい、hg19 or GRCh38"       )
p <- add_argument(p, "--chr",              default='chr12',             type="character", help="chr-start-endの領域をプロットする"                                       )
p <- add_argument(p, "--start",            default=120291763,           type="numeric",   help="chr-start-endの領域をプロットする"                                       )
p <- add_argument(p, "--end",              default=120291903,           type="numeric",   help="chr-start-endの領域をプロットする"                                       )
p <- add_argument(p, "--domain_file",      default='rnu4-2_domain.txt', type="character", help="name,start,endの名前の3列でドメイン領域を定義したタブ区切りファイル"     )
p <- add_argument(p, "--out_png",          default='output.png',        type="character", help=""                                                                        )
argv <- parse_args(p)


CHR = argv$chr
STA = argv$start
END = argv$end


# 領域を指定してそこだけVCFを読み込む
region_param <- ScanVcfParam(which=GRanges(CHR, IRanges(STA, END)))
vcf <- readVcf(argv$clinvar_vcf, argv$reference_genome, param=region_param)


# 必要な情報をVCFから取ってデータフレームにする
# CLNSIG, GENEINFO, CLNHGVSなど
# で病原性で色分けする
chrom <- seqnames(rowRanges(vcf))
pos <- start(rowRanges(vcf))
info_df <- info(vcf)
clnsig <- info_df$CLNSIG
geneinfo <- info_df$GENEINFO
clnhgvs <- info_df$CLNHGVS

df <- tibble(
  CHROM = as.character(chrom),
  POS = pos,
  CLNSIG = unlist(clnsig),
  GENEINFO = geneinfo,
  CLNHGVS = unlist(clnhgvs)
)

df$col <- case_when(
  df$CLNSIG == "Benign" ~ "blue",
  df$CLNSIG == "Likely_benign" ~ "blue",
  df$CLNSIG == "Uncertain_significance" ~ "yellow",
  df$CLNSIG == "Likely_pathogenic" ~ "red",
  df$CLNSIG == "Pathogenic/Likely_pathogenic" ~ "red",
  df$CLNSIG == "Pathogenic" ~ "red",
  TRUE ~ "gray"
)


# で、そのdfをlolliplotに読み込ませる用に変える
sample.gr <- GRanges( CHR, IRanges(df$POS, width=1, names=df$CLNHGVS))
sample.gr$color <- df$col


# ドメイン領域（とその色と高さ）を定義する
domains <- read.table(argv$domain_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

features <- GRanges(
  CHR,
  IRanges(
    start = domains$start,
    end = domains$end,
    names = domains$name
  )
)

features$fill = rep("#98CE31", 5)
features$height <- 0.05


# lolliplotで描いて保存する
xaxis <- c(STA, END)
ranges <- GRanges(CHR, IRanges(STA, END))

legend <- c(
  "blue",
  "yellow",
  "red"
)

names(legend) <- c(
  "Benign or Likely_benign",
  "Uncertain_significance",
  "Pathogenic or Likely_pathogenic"
)

png(argv$out_png, units='cm', width=20, height=10, res=600)
lolliplot(sample.gr, features, xaxis=xaxis, ranges=ranges, legend=legend, cex=.4)
dev.off()



#dict_keys(['obj_type', 'accession', 'accession_version', 'title', 'variation_set', 'supporting_submissions', 'germline_classification', 'clinical_impact_classification', 'oncogenicity_classification', 'record_status', 'gene_sort', 'chr_sort', 'location_sort', 'variation_set_name', 'variation_set_id', 'genes', 'molecular_consequence_list', 'protein_change', 'fda_recognized_database'])
#
#
#
#from Bio import Entrez
#Entrez.email = "your.email@example.com"
#
## 1. esearch
#search_handle = Entrez.esearch(db="clinvar", term="TP53[gene]", retmax=10)
#search_results = Entrez.read(search_handle)
#ids = search_results["IdList"]
#
## 2. esummary（validate=False に変更）
#summary_handle = Entrez.esummary(db="clinvar", id=",".join(ids), retmode="xml")
#summary_records = Entrez.read(summary_handle, validate=False)
#
## 3. 情報の抽出
#for record in summary_records['DocumentSummarySet']['DocumentSummary']:
#    uid = record.attributes['uid']
#    accession = record.get('accession', 'N/A')
#    variation_name = record.get('variation_name', 'N/A')
#    print(f"{uid}\t{accession}\t{variation_name}")
#
#
#
#import pandas as pd
#
#records = summary_records['DocumentSummarySet']['DocumentSummary']
#
#def flatten_record(record):
#    # recordはBio.Entrezのパーサー結果のXMLノードなのでdict-likeだが、
#    # ネストしたリストや辞書があることが多い
#    flat = {}
#    # まずはトップレベルのキーと値を追加（値が文字列か数字の場合のみ）
#    for key, value in record.items():
#        # もし値がリストや辞書なら一旦スキップするか文字列化する
#        if isinstance(value, (str, int, float)):
#            flat[key] = value
#        elif value is None:
#            flat[key] = None
#        else:
#            # 複雑なオブジェクトは文字列化で一旦保存
#            flat[key] = str(value)
#    #    
#    # 例えば 'genes' はリストなので最初の遺伝子名だけ抜き出し（例）
#    genes = record.get('genes')
#    if genes and isinstance(genes, list):
#        flat['genes_first'] = genes[0] if len(genes) > 0 else None
#    else:
#        flat['genes_first'] = None
#    #    
#    # 'protein_change' も複数あるかもしれないので一旦文字列化して格納
#    protein_change = record.get('protein_change')
#    if protein_change:
#        if isinstance(protein_change, list):
#            flat['protein_change'] = ";".join(str(pc) for pc in protein_change)
#        else:
#            flat['protein_change'] = str(protein_change)
#    else:
#        flat['protein_change'] = None
#    #    
#    # 'variation_set' はリストや複雑なネストがあるので文字列化する例
#    variation_set = record.get('variation_set')
#    if variation_set:
#        flat['variation_set'] = str(variation_set)
#    else:
#        flat['variation_set'] = None
#    #    
#    return flat
#
## すべてのレコードを展開してリストに
#flat_records = [flatten_record(rec) for rec in records]
#
## DataFrameに変換
#df = pd.DataFrame(flat_records)
#
#print(df.head())
#
#
