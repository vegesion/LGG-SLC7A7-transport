# =============================================================================
# 11_depth_sensitivity.R                 [Supplementary — 리뷰어 방어이자 차별점]
# 세포 단위 +/- 이분의 depth 교란 정량화: depth 비 · glmer OR · 매칭 음성대조
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(lme4); library(AUCell) })

mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
sets <- list(transport = build_transporter_geneset(universe = rownames(mg)),
             oxphos = get_msig_genes("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
             tnfa   = get_msig_genes("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
auc <- score_auc(mg, sets)
for (nm in colnames(auc)) mg[[paste0("AUC_", nm)]] <- auc[, nm]

res <- dplyr::bind_rows(lapply(colnames(auc), function(nm) {
  s <- depth_sensitivity(mg, module_col = paste0("AUC_", nm))
  data.frame(module = nm, depth_ratio = s$depth_ratio,
             glmer_OR = if (is.null(s$glmer)) NA else s$glmer$OR,
             glmer_p  = if (is.null(s$glmer)) NA else s$glmer$p,
             observed_delta = s$observed_delta,
             null_median = stats::median(s$null_delta, na.rm = TRUE),
             empirical_p = s$empirical_p)
}))
print(res); write_result(res, "Supp_depth_sensitivity.csv")
# 해석: depth_ratio 가 1에서 멀수록 +/- 가 depth 로 구분됨. glmer_OR 이 유의하면
#       "SLC7A7 검출 = depth" 교란의 직접 증거. empirical_p 가 크면 그 모듈 차이는
#       SLC7A7 특이적이라 보기 어렵다(음성대조 분포 안).
message("완료: depth 민감도")
