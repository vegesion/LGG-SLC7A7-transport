# =============================================================================
# 01_tcga_download_prepare.R
# TCGA-LGG / TCGA-GBM STAR-counts 를 "한 번만" 다운로드 → 로컬 저장 → logCPM 정규화
# 입력: GDC (온라인, 최초 1회)   출력: data/processed/{LGG,GBM}_* , data/raw/GDCdata/
# 원본: TCGA.R, TCGA_GBM.R, egdeR.R
# -----------------------------------------------------------------------------
# ★ 캐시 방식: get_tcga_data() 는 이미 저장된 .rds 가 있으면 그걸 읽고,
#   없을 때만 GDC에서 받습니다. 즉 이 스크립트를 다시 돌려도 재다운로드하지 않습니다.
#   원본 파일은 data/raw/GDCdata/ 에 영구 보관됩니다(수백 개 폴더가 생기는 건 정상).
#   강제로 다시 받으려면 각 호출에 force = TRUE.
#   다운로드는 실패 시 자동 재시도(최대 3회)하며, 다운로드 일시는
#   data/processed/{tag}_download_info.txt 와 rds의 out$download_info 에 기록됩니다.
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment)
  library(AnnotationDbi); library(org.Hs.eg.db); library(edgeR)
})

# 큰 파일 다운로드가 타임아웃에 잘리지 않도록 (재시도 루프 예방)
options(timeout = 10000)

for (proj in c("TCGA-LGG", "TCGA-GBM")) {
  tag <- sub("TCGA-", "", proj)

  # 1) 다운로드/준비 (캐시 우선) → list(counts = SYMBOL count, clinical)
  out <- get_tcga_data(proj)                     # 이미 있으면 재다운로드 안 함

  # 2) 원시 count(txt) 백업 (없을 때만)
  raw_txt <- file.path(DIR_DATA_PROC, paste0(tag, "_TCGA_rawcounts.txt"))
  if (!file.exists(raw_txt)) {
    write.table(out$counts, raw_txt, sep = "\t", quote = FALSE, fileEncoding = "UTF-8")
  }

  # 3) logCPM 정규화 (edgeR::cpm, log=TRUE, prior.count=1) → norm csv (없을 때만)
  norm_csv <- file.path(DIR_DATA_PROC, paste0(tag, "_TCGA_norm.csv"))
  if (!file.exists(norm_csv)) {
    norm <- normalize_logcpm(out$counts, prior_count = 1)
    write.csv(norm, norm_csv, row.names = FALSE)
    message("[saved] ", norm_csv)
  } else {
    message("[cache] norm csv 이미 있음 → ", norm_csv)
  }
}
message("완료: TCGA 준비. 재실행 시 로컬 저장본을 사용하므로 재다운로드하지 않습니다.")
