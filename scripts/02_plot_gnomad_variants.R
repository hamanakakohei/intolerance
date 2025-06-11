#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(patchwork)
library(VariantAnnotation)
library(argparser)

source("scripts/intolerance_utils.R")


###### 方針、やること######
# gencodeのhuman gtfもしくはgffをAPI？REST？か何かでアクセスして、
# 指定したトランスクリプトの配列を得たい。もちろん複数のエクソンを一つにつなげたい。
# だが、これは気軽に少数遺伝子を見たい時の話で、データベースを作るならちゃんと元ファイルをダウンロードすべき
# 
# トランスクリプト配列
# mutation
# - snv, del, insまとめてロリポップでok
# alt allele n (snv)
# - 閾値を変えて棒グラフ（max Y=3）を何パターンか作る
# -- ACが1以上
# -- ACが2以上
# -- ACが3以上
# alt allele n (snv + del[範囲を考慮])
# - 閾値を変えて棒グラフ（max Y=4）を何パターンか作る
# -- ACが1以上
# -- ACが2以上
# -- ACが3以上
# alt allele n (ins)
# - 閾値を変えて棒グラフ（max Y=無限、insの中身次第）を何パターンか作る
# -- ACが1以上
# -- ACが2以上
# -- ACが3以上
# agg. minor AC (snv)
# - 1行ヒートマップかrectangleで作る
# - 0, 1, 2, 3, >3の5色くらいで十分
# agg. minor AC (snv + del[範囲を考慮])
# - 1行ヒートマップかrectangleで作る
# - 0, 1, 2, 3, >3の5色くらいで十分
# each minor AC (snv or snv + del)
# - 4行ヒートマップ（del含めるなら5行）
# - refと同じ塩基は灰色
# - alt塩基はACで色分けするが、0, 1, 2, 3, >3の5色くらいで十分
# - インプットはnormalized VCF
# - スクリプトはどこかから取ってきたい
# each minor AC (ins)
# - 積み上げ棒グラフ
# - 各マスの中に挿入長を書く、2桁以上はアスキー文字？
# - ACで色分けするが、0, 1, 2, 3, >3の5色くらいで十分
#############################


p <- arg_parser("")
p <- add_argument(p, "--chr",     default='chr12',      help="chr-start-endの領域をプロットする",  type="character")
p <- add_argument(p, "--start",   default=120291763,    help="chr-start-endの領域をプロットする",  type="numeric")
p <- add_argument(p, "--end",     default=120291903,    help="chr-start-endの領域をプロットする",  type="numeric")
p <- add_argument(p, "--ref_ver", default='GRCh38',     help="VariantAnnotationパッケージの挙動に影響するらしい",  type="character")
p <- add_argument(p, "--vcf",     default='output.vcf', help="get_gene-level_gnomad_info.pyで得たvcf、normalizedしなくても大丈夫なはず、、",  type="character")
p <- add_argument(p, "--min_ac",  default=1,            help="alt allele種類数ヒストグラムを描く時のalt allele ACの閾値",  type="numeric")
p <- add_argument(p, "--out",     default='output.png', help="",  type="character")
argv <- parse_args(p)

CHR = argv$chr
STA = argv$start
END = argv$end


# 領域を絞る
vcf <- readVcf(argv$vcf, argv$ref_ver)
target_gr <- GRanges( CHR, IRanges( STA, END ))
vcf <- vcf[overlapsAny(rowRanges(vcf), target_gr)]


# SNV のみに絞る
vcf <- expand(vcf)
is_snv <- width(ref(vcf)) == 1 & width(alt(vcf)) == 1
vcf <- vcf[is_snv]


# 全ポジ x 全塩基のACのデータフレームを作る
all_pos <- STA:END
alt_alleles <- c("A", "C", "G", "T")

all_combos_df <- expand.grid(
    pos = all_pos,
    alt = alt_alleles,
    stringsAsFactors = FALSE
  ) %>% 
  as_tibble

all_pos_df <- tibble(pos = all_pos)

ac_df <- tibble(
  pos = start(rowRanges(vcf)),
  alt = as.character(alt(vcf)),
  ac = info(vcf)$AC
  #ac = unlist( info(vcf)$AC )
)

ac_heatmap_df <- left_join(all_combos_df, ac_df, by = c("pos", "alt")) %>%
  mutate(ac = ifelse(is.na(ac), 0, ac))

ac_heatmap_df <- ac_heatmap_df %>%
  mutate(value = case_when(
    is.na(ac) ~ "0",
    ac == 0 ~ "0",
    ac == 1 ~ "1",
    ac == 2 ~ "2",
    ac >= 3 ~ "3+"
  ))


# 全ポジのアグリゲートACのデータフレームを作る
ac_heatmap_df_agg <- ac_heatmap_df %>%
  group_by(pos) %>%
  summarise(ac = sum(ac, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(alt='agg')

ac_heatmap_df_agg <- ac_heatmap_df_agg %>%
  mutate(value = case_when(
    is.na(ac) ~ "0",
    ac == 0 ~ "0",
    ac == 1 ~ "1",
    ac == 2 ~ "2",
    ac >= 3 ~ "3+"
  ))


# 全ポジのalt allele種類数のデータフレームを作る
alt_n_df <- ac_df %>%
  group_by(pos) %>%
  summarise(distinct_alt_n = n_distinct(alt), .groups = "drop") %>%
  right_join(all_pos_df) %>%
  mutate(distinct_alt_n = ifelse(is.na(distinct_alt_n), 0, distinct_alt_n))


# 図を描いて結合して保存
p_bar1_1 <- make_bar(alt_n_df, "Distinct alt allele (SNV, minor AC > 0)", "gray30")
#p_bar1_2 <- make_bar("Distinct alt allele (SNV, minor AC > 1)", "gray30")
#p_bar1_3 <- make_bar("Distinct alt allele (SNV, minor AC > 2)", "gray30")
#p_bar2_1 <- make_bar("Distinct alt allele (SNV & DEL, minor AC > 0)", "gray50")
#p_bar2_2 <- make_bar("Distinct alt allele (SNV & DEL, minor AC > 1)", "gray50")
#p_bar2_3 <- make_bar("Distinct alt allele (SNV & DEL, minor AC > 2)", "gray50")
#p_bar3_1 <- make_bar("Distinct alt allele (INS, minor AC > 0)", "gray70")
#p_bar3_2 <- make_bar("Distinct alt allele (INS, minor AC > 1)", "gray70")
#p_bar3_3 <- make_bar("Distinct alt allele (INS, minor AC > 2)", "gray70")

p_hm4 <- make_heatmap(ac_heatmap_df, "Minor AC (SNV)")
p_hm1_1 <- make_heatmap(ac_heatmap_df_agg, "Aggregated minor AC (SNV)", show_legend=TRUE)
#p_hm5 <- make_heatmap(hm5, "Minor AC (SNV & DEL)", show_legend=FALSE)
#p_hm1_2 <- make_heatmap(hm1_2, "Aggregated minor AC (SNV & DEL)", show_legend=TRUE)

#p <- p_bar1_1 / p_bar1_2 / p_bar1_3 / p_bar2_1 / p_bar2_2 / p_bar2_3 / p_bar3_1 / p_bar3_2 / p_bar3_3 / p_hm1_1 / p_hm1_2 / p_hm4 / p_hm5 / p_bottom +
#plot_layout(heights = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 0.5, 2, 2, 1))  # 各段の高さ比
p <- p_bar1_1 / p_hm1_1 / p_hm4 +
  plot_layout(heights = c(1, 0.5, 2))  # 各段の高さ比

ggsave( argv$out, plot = p, width = 20, height = 10, dpi = 600, units = "cm")



#### 1. INS用のカラフル棒グラフ（上下反転）
#
#df <- data.frame(
#  bin = rep(c("Bin1", "Bin2", "Bin3", "Bin4"), times = c(4, 7, 1, 6)),
#  y = c(1:4, 1:7, 1, 1:6),
#  label = c(strsplit("ABCD", "")[[1]],
#            strsplit("EFGHIJK", "")[[1]],
#            "L",
#            strsplit("MNOPQR", "")[[1]]),
#  fill = c("red", "orange", "red", "orange",           
#           "blue", "cyan", "blue", "cyan", "blue", "cyan", "blue",  
#           "green",                                     
#           "purple", "violet", "purple", "violet", "purple", "violet")  
#)
#
#p_bottom <- ggplot() +
#  geom_col(data = df, aes(x = bin, y = 1, fill = fill),
#           position = position_stack(vjust = 0), width = 1) +
#  geom_text(data = df, aes(x = bin, y = y - 0.5, label = label),
#            color = "white", size = 5, fontface = "bold") +
#  scale_fill_identity() +
#  #scale_y_reverse() +  # ここで上下反転
#  theme_minimal() +
#  labs(title = NULL, x = NULL, y = NULL) +
#  theme(axis.text.x = element_blank(),
#        axis.title.y = element_blank(),
#        plot.margin = margin(5, 5, 5, 5))
