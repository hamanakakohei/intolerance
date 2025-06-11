#!/usr/bin/env Rscript

library(VariantAnnotation)
library(trackViewer)
library(dplyr)
library(argparser)

p <- arg_parser("")
p <- add_argument(p, "--clinvar_vcf",      default='clinvar.vcf.gz',    type="character", help="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gzとtbi" )
p <- add_argument(p, "--reference_genome", default='hg38',              type="character", help="VariantAnnotationパッケージの挙動に影響するらしい"       		 )
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
