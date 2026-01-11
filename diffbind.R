BiocManager::install("DiffBind")
library(DiffBind)
library(readxl)
library(BiocParallel)
# Windows下仅能开启“伪并行”（进程模拟，提速有限）
register(SnowParam(workers = 4))  # 4是线程数
setwd("C:/Users/19892/Desktop/ATAC")
group<-read_excel("ATAC_group.xlsx")
group<-group[,3:4]
setwd("C:/Users/19892/Desktop/ATAC/PRJNA439280")
# === 处理 GSE210285 ===
peaks_dir <- file.path("peaks")
bam_dir   <- file.path("bam_filtered")
# 获取 peak 文件和样本名
peak_files <- list.files(peaks_dir, pattern = "_peaks.narrowPeak$", full.names = TRUE)
sample_names <- gsub("_peaks.narrowPeak$", "", basename(peak_files))
# 构建 sample sheet（必须包含 SampleID, Condition, bamReads, Peaks）
samples<- data.frame(
  SampleID = sample_names,
  bamReads = file.path(bam_dir, paste0(sample_names, ".clean.bam")),
  Peaks = peak_files,
  stringsAsFactors = FALSE
)
samples<-merge(samples,group,by.x="SampleID",by.y="SRR_id")

samples$Factor<-as.factor(samples$Factor)
# 检查是否匹配
print(samples)
db <- dba(sampleSheet=samples, scoreCol = 5)
# 计数（在 consensus peak 区域内统计 reads）
db <- dba.count(db,bParallel =T)
# 设置对比：比如 Treatment vs Control
db_contrast <- dba.contrast(db, contrast =c("Factor","RASOIS_ATACSEQ_T0","RASOIS_ATACSEQ_144H"))
# 进行差异分析
db_analyzed <- dba.analyze(db_contrast,method=DBA_DESEQ2)
#  summary of results
dba.show(db_analyzed, bContrasts=T)
# 获取所有结果（不设阈值）
res <- dba.report(db_analyzed)  # th=1 表示返回所有 peaks
res<-as.data.frame(res)
write.table(res,"RASOIS_ATACSEQ_144H_vs_RASOIS_ATACSEQ_T0_peaks_diff_results_fdr0.05_Deseq2.txt")
