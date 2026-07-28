# =============================================================================
# R/setup.R  ---  세션 초기화 / session bootstrap
# -----------------------------------------------------------------------------
# 모든 analysis/ 스크립트 첫 줄에서  source(here::here("R", "setup.R"))  로 호출.
# config + 모든 함수 모듈 + 공통 라이브러리를 한 번에 로드합니다.
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

# 설정 로드
source(here::here("config", "config.R"))

# 공통 라이브러리 (분석 전반 공유). 도메인 특화 패키지는 각 스크립트에서 로드.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

# 함수 모듈 로드 (순서 무관, 상호 참조 없음)
for (f in c("utils.R", "tcga_data.R", "survival.R", "gsea.R", "scrna.R",
             "gene_sets.R", "cohorts.R", "signature.R", "validation.R",
             "deconvolution.R", "sc_donor.R", "spatial.R", "mr.R")) {
  source(here::here("R", f))
}

message("[setup] 프로젝트 로드 완료 | 대상 유전자: ", GENE_OF_INTEREST,
        " | 루트: ", here::here())
