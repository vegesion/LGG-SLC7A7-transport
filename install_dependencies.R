# =============================================================================
# install_dependencies.R  ---  최초 1회 실행 / one-time dependency install
# -----------------------------------------------------------------------------
# ※ "Error in readRDS(dest): 알 수 없는 입력 포맷" 은 저장소 메타데이터/패키지
#   다운로드가 깨졌을 때 나는 오류입니다. 아래 설정(libcurl + 넉넉한 timeout +
#   올바른 Bioconductor 저장소)이 그 문제를 예방합니다.
#   그래도 나면: (1) 새 R 세션에서 재시도  (2) 사내 프록시/VPN 이면 잠시 끄기
#   (3) 다른 CRAN 미러 사용.
# =============================================================================

# ── 0. 안정적 다운로드 환경 ──────────────────────────────────────────────────
options(
  repos                = c(CRAN = "https://cloud.r-project.org"),
  download.file.method = "libcurl",   # Windows 다운로드 깨짐(readRDS 오류) 예방
  timeout              = 600          # 큰 메타/패키지 파일 대비 (기본 60초는 짧음)
)

# ── 1. BiocManager + Bioconductor 3.20 저장소 확립 ──────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 이 한 줄이 Bioconductor 저장소를 repos 에 제대로 등록합니다(경고 해소).
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

# ── 2. 한 개 실패해도 나머지는 계속 설치하는 헬퍼 ───────────────────────────
install_safe <- function(pkgs, installer) {
  need <- setdiff(pkgs, rownames(installed.packages()))
  if (!length(need)) { message("모두 설치됨."); return(invisible()) }
  failed <- character(0)
  for (p in need) {
    message(">>> 설치: ", p)
    ok <- tryCatch({ installer(p); TRUE },
                   error = function(e) { message("  실패: ", p, " — ", conditionMessage(e)); FALSE })
    if (!ok) failed <- c(failed, p)
  }
  if (length(failed)) message("\n[미설치] ", paste(failed, collapse = ", "),
                              "\n  → 새 세션에서 개별 재시도하세요.")
  invisible(failed)
}

# ── 3. CRAN 패키지 ───────────────────────────────────────────────────────────
cran <- c("here", "tidyverse", "survival", "survminer", "coxphf", "devEMF",
          "pbapply", "ggpubr", "ggrepel", "rstatix", "viridis", "patchwork",
          "jsonlite", "remotes", "devtools", "anndata", "reticulate", "readr", "hdf5r",
          "glmnet", "timeROC", "lme4", "lmerTest", "pheatmap")
install_safe(cran, function(p) install.packages(p))

# ── 4. Bioconductor 패키지 ───────────────────────────────────────────────────
bioc <- c("TCGAbiolinks", "SummarizedExperiment", "AnnotationDbi", "org.Hs.eg.db",
          "clusterProfiler", "msigdbr", "enrichplot", "fgsea", "limma", "edgeR",
          "DESeq2", "Seurat", "SingleCellExperiment", "slingshot", "biomaRt", "UCell")
install_safe(bioc, function(p) BiocManager::install(p, ask = FALSE, update = FALSE))

# ── 5. GitHub 전용 패키지 ────────────────────────────────────────────────────
gh <- c(BPCells     = "bnprks/BPCells/r",
        xCell       = "dviraran/xCell",
        CellChat    = "jinworks/CellChat",
        presto      = "immunogenomics/presto",
        monocle3    = "cole-trapnell-lab/monocle3",
        MRcML       = "xue-hr/MRcML",
        TwoSampleMR = "MRCIEU/TwoSampleMR",
        MRPRESSO    = "rondolab/MR-PRESSO")
install_safe(names(gh), function(p) remotes::install_github(gh[[p]], upgrade = "never"))


# # -----estimate 설치------ ------------------------------------------------


install.packages("estimate", repos="http://R-Forge.R-project.org", dependencies=TRUE)

message("\n의존성 설치 절차 완료. 미설치 항목이 있으면 위 로그를 확인하세요.")
