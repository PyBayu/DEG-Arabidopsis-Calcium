#################################################################
# Analisis Ekspresi Gen Ca
# Dataset: GSE278279 (Normal Ca vs Low Ca)
# Platform: Microarray (Affymetric GPL198)
# Tujuan : Mengidentifikasi Differentially Expressed Genes (DEG)
#################################################################

#################################################################
# PART A. PENGANTAR KONSEP
#################################################################

# Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen
# antara dua kondisi biologis (Normal Ca dan Low Ca).
# Pada skrip ini akan menggunakan pendekatan statistic limma (Linear Models for Microarray Data), yang merupakan standar emas untuk data microarray.

################################
# PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & lOAD PACKAGE)
################################

# ---- Apa itu package? ----
# Package adalah kumpulan fungsi siap pakai di R.
# Bioinformatika di R sangat bergantung pada package dari CRAN dari Bioconductor.

# 1. Instal BiocManager (manajer paker Bioconductor)
# IF adalah struktur logika: " jika kondisi terpenuhi, lakukan aksi"
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}

# 2. Instal paket Bioconductor
# GEOquery: mengambil data dari database GEO
# limma : analisis statistic ekspresi gen
BiocManager::install(c(
  "GEOquery",
  "limma"
), ask = FALSE, update = FALSE)
# Install annotation package sesuai platform
# GPL198 =  Affymetrix Arabidopsis ATH1 Genome Array
BiocManager::install("ath1121501.db",
                     ask = FALSE, update = FALSE)

# 3. Install paket CRAN untuk visualisasi dan manipulasi data
# pheatmap: heatmap ekspresi gen
# ggplot2: grafik (volcano plot)
#dplyr : manipulasi tabel data
# install.packages(c(pheatmap", "ggplot2", "dplyr"))
# umap: grafik (plot UMAP)
# if (!requireNamespace("umap", quietly = TRUE)){install.packages("umap")}

# 4. Memanggil library
# library() digunaka agar fungsi di dalam package bisa digunakan
library (GEOquery)
library(Biobase)
library(BiocGenerics)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(ath1121501.db)
library(umap)

#################################################################
# PART C. PENGAMBILAN DATA DARI GEO
#################################################################

# GEO (Gene Expression Omnibus) adala database public milik NCBI
# getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
# GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
# AnnotGPL = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE278279",
               GSEMatrix = TRUE,
               AnnotGPL = TRUE)[[1]]

# ExpressionSet berisi:
# -exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

################################
# PART D. PRE-PROCESSING DATA EKSPRESI
################################

# exprs(): mengambil matriks ekspresi gen
# Baris = probe/gen
# Kolom = sampel
ex <- exprs(gset)
View(ex)
# ---- Mengapa perlu log2 transfromasi? ----
# Data microarray mentah memiliki rentang nilai sangat besar.
# Log2 digunakan untuk:
# 1. Menstabilkan varians
# 2. Mendekati asumsi model linear
# 3. Memudahkan interpretasi log fold change

# quantile(): menghitung nilai kuantil (persentil)
# as.numeric(): mengubah hasil quantile ( yang berupa named vector) menjadi vector bumerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

# LogTransform adalah variable logika (TRUE/FALSE)
# Operator logika:
# > : lebih besar dari
# || : OR (atau)
# && : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx [2] > 0)

# IF statement:
# Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  
  # Transformasi log2
  ex <- log2(ex)
}

#################################################################
# PART E. DEFINISI KELOMPOK SAMPEL
#################################################################

# pData(): metadata sampel
# source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]

# make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names (group_info)

#factor():
# Mengubah data kategorik menjadi factor
# Faktor sangat penting untuk analisis statistic di R
gset$group <- factor(groups)

# levels(): melihat kategori unik dalam factor
nama_grup <-levels(gset$group)
print(nama_grup)

#Khusus
#################################################################
# Ubah 4 grup jadi 2 grup
#################################################################
# Ambil informasi kondisi sampel
group_info <- pData(gset)[["source_name_ch1"]]

# Kelompokkan berdasarkan kadar Ca
groups <- ifelse(grepl("0.3", group_info),
                 "LowCa",
                 "NormalCa")

# Ubah menjadi factor
gset$group <- factor(groups)

# Cek hasil
levels(gset$group)

## Asli
#################################################################
# PART F. DESIGN MATRIX (KERANGKA STATISTIK)
#################################################################

# model.matrix():
# Membuat matriks desain untuk model linear
# ~0 berarti TANPA intercept (best practice limma)

design <- model.matrix(~0 + gset$group)

# colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels (gset$group)

# Menentukan perbandingan biologis
grup_normal_ca <- nama_grup[1]
grup_low_ca <- nama_grup[2]

contrast_formula <- paste(grup_normal_ca, "-", grup_normal_ca)
print(paste("Kontras yang dianalisis:", contrast_formula))

# Khusus Ca
#################################################################
# PART F. DESIGN MATRIX
#################################################################

# Membuat design matrix
design <- model.matrix(~0 + gset$group)

# Memberi nama kolom
colnames(design) <- levels(gset$group)

# Melihat nama grup
nama_grup <- levels(gset$group)

# Menentukan grup
grup_low_ca <- nama_grup[1]
grup_normal_ca <- nama_grup[2]

# Membuat kontras
contrast_formula <- paste(grup_low_ca, "-", grup_normal_ca)

print(paste("Kontras yang dianalisis:",
            contrast_formula))
View(design)

#################################################################
# PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
#################################################################

# lmFit():
# Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

# makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

#eBayes():
# Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

# topTable()
# Mengambil hasil akshir DEG
# adjust = "fdr" -> koreksi multiple testing
# p.value = 0,01 -> gen sangat signifikan
topTableResults <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf, p.value =0.01)

head(topTableResults)

####################################
# PART H. ANOTASI NAMA GEN
####################################

# Penting: # Pada data microarray Affymetrix, unit analisis awal adalah PROBE,
# bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan database resmi Bioconductor.

# Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

# Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(ath1121501.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")

# Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(topTableResults, gene_annotation, by = "PROBEID", all.x = TRUE)

# Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])
View(topTableResults)

#################################################################
# PART I.1. BOXPLOT DISTRIBUSI NILAI EKSPRESI
#################################################################

# Boxplot digunakan untuk:
# - Mengecek distribusi nilai ekspresi antar sampel
# - Melihat apakah ada batch effect
# - Mengevaluasi apakah normalisasi/log-transform wajar

# Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)
boxplot(ex,
        col = group_colors,
        las = 2,
        outline = FALSE,
        main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
        ylab = "Expression Value (log2)",
        cex.axis = 0.75)

legend("topright",
       legend = levels(gset$group),
       fill = unique(group_colors),
       cex = 0.8)

#################################################################
# PART I.2. DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
#################################################################

# Density plot menunjukkan sebara global nilai ekspresi gen
# Digunakan untuk:
# - Mengecek efek log-transform
# - Membandingkan distribusi antar grup

# Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long,
       aes(x = Expression,
           color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (log2)",
    y = "Density"
  )

#################################################################
# PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)
#################################################################

# UMAP digunakan untuk:
# - Mereduksi ribuan menjadi 2 dimensi
# - Melihat pemisahan sampel secara global
# - Alternatif PCA (lebih sensitive ke struktur local)

# Transpose matriks ekspresi:
# UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

# Jalankan UMAP
umap_result <- umap(umap_input)

# Versi lain
n_sample <- nrow(umap_input)

umap_result <- umap(
  umap_input,
  n_neighbors = min(5, n_sample - 1)
)

# Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Group = gset$group
)

# Plot UMAP
ggplot(umap_df,
       aes(x = UMAP1,
           y = UMAP2,
           color = Group)) +
  geom_point(size = 3,
             alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#####################################
# PART J.1. VISUALISASI VOLCANO PLOT
#####################################

# Volcano plot menggabungkan:
# - Log fold change (efek biologis)
# - Signifikansi statistic

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

# Klasifikasi status gen
volcano_data$status <- "NO"

volcano_data$status[
  volcano_data$logFC > 1 &
    volcano_data$adj.P.Val < 0.01
] <- "UP"

volcano_data$status[
  volcano_data$logFC < -1 &
    volcano_data$adj.P.Val < 0.01
] <- "DOWN"

# Visualisasi
ggplot(volcano_data,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = status)) +
  
  geom_point(alpha = 0.6) +
  
  scale_color_manual(values = c(
    "DOWN" = "blue",
    "NO" = "grey",
    "UP" = "red"
  )) +
  
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  
  theme_minimal() +
  
  ggtitle("Volcano Plot DEG Respon Tanaman Terhadap Kalsium")


#######################################
# PART J.2. VISUALISASI HEATMAP
#####################################

# Heatmap digunakan untuk melihat pola ekspresi gen antar sampel berdasarkan gen-gen paling signifikan

# Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

# Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

# Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "", 
  top50$PROBEID, # jika SYMBOL kosong -> probe ID
  top50$SYMBOL # jika ada -> gene symbol
)

rownames(mat_heatmap) <- gene_label

# Pembersihan data (WAJIB agar tidak error hclust)
# Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

# Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)

mat_heatmap <- mat_heatmap[
  gene_variance > 0,
]


# Versi lain
# Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

# Paksa jadi matriks numerik
mat_heatmap <- as.matrix(mat_heatmap)
mode(mat_heatmap) <- "numeric"

# Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
class(mat_heatmap)
var
rm(var)
var
# 5. Anotasi kolom (kelompo sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

# Visualisasi heatmap (final, tidak crowded)
pheatmap(
  mat_heatmap,
  
  scale = "row", # z-score per gen
  
  annotation_col = annotation_col,
  
  show_colnames = FALSE, #nama sampel dimatikan
  
  show_rownames = TRUE,
  
  fontsize_row = 7,
  
  clustering_distance_rows = "euclidean",
  
  clustering_distance_cols = "euclidean",
  
  clustering_method = "complete",
  
  main = "Top 50 Differentially Expressed Genes"
)

##########################
# PART  K. MENYIMPAN HASIL
############################

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE278279_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")

##########################
# PART K. MENYIMPAN HASIL
##########################

# Tentukan folder tujuan
folder_tujuan <- "D:/Baru"

# Simpan hasil DEG ke CSV
write.csv(
  topTableResults,
  file = paste0(folder_tujuan, "/Hasil_GSE278279_DEG.csv"),
  row.names = FALSE
)

message("Analisis selesai. File hasil telah disimpan di D:/Baru")

# Penyimpan up down di excel
library(openxlsx)
# Klasifikasi DEG
topTableResults$DEG_status <- "NO"

topTableResults$DEG_status[
  topTableResults$logFC > 1 &
    topTableResults$adj.P.Val < 0.01
] <- "UP"

topTableResults$DEG_status[
  topTableResults$logFC < -1 &
    topTableResults$adj.P.Val < 0.01
] <- "DOWN"

# Pisahkan Up dan Dwon
deg_up <- topTableResults %>%
  filter(DEG_status == "UP")

deg_down <- topTableResults %>%
  filter(DEG_status == "DOWN")

# Buat ringkasan jumlah gen
summary_deg <- data.frame(
  Category = c("Total_genes", "UP_genes", "DOWN_genes"),
  Count = c(
    nrow(topTableResults),
    nrow(deg_up),
    nrow(deg_down)
  )
)

# Simpan ke excel
# Folder tujuan
folder_tujuan <- "D:/Baru"

# Nama file
file_excel <- paste0(folder_tujuan, "/Hasil_DEG_GSE278279.xlsx")

# Buat workbook
wb <- createWorkbook()

# Tambah sheet
addWorksheet(wb, "All_DEG")
addWorksheet(wb, "UP_genes")
addWorksheet(wb, "DOWN_genes")
addWorksheet(wb, "Summary")

# Isi data
writeData(wb, "All_DEG", topTableResults)
writeData(wb, "UP_genes", deg_up)
writeData(wb, "DOWN_genes", deg_down)
writeData(wb, "Summary", summary_deg)

# Simpan file
saveWorkbook(wb, file_excel, overwrite = TRUE)

message("File Excel berhasil disimpan di D:/Baru")
