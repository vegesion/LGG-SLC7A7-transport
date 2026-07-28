# =============================================================================
# 00_donor_pseudobulk_icc.R   ★ 관문(gate)
# GBmap 배치 구조 점검: SLC7A7 발현 분산이 donor 에서 오는가 dataset 에서 오는가?
# 입력: data/processed/microglia.rds (06)   출력: results/ICC_*.csv
# ※ lme4::lmer 는 p-value 를 주지 않음 → lmerTest 가 있으면 사용(Satterthwaite).
# =============================================================================
source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(lme4) })
if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)

obj <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))
ds_col <- find_batch_col(obj)
if (is.na(ds_col)) stop("배치 컬럼을 찾지 못했습니다: ", paste(colnames(obj@meta.data), collapse = ", "))
message("배치 컬럼: ", ds_col)

df <- data.frame(expr = FetchData(obj, vars = GENE_OF_INTEREST)[, 1],
                 depth = obj[[DEPTH_COVARIATE]][, 1],
                 donor = obj$donor_id[, 1],
                 dataset = as.character(obj[[ds_col]][, 1]))
agg <- df %>% dplyr::group_by(dataset, donor) %>%
  dplyr::summarise(expr = mean(expr), depth = mean(depth), n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= SC_DONOR_MIN_CELLS)
message("donor ", nrow(agg), " / dataset ", dplyr::n_distinct(agg$dataset))

m  <- lmer(expr ~ scale(depth) + (1 | dataset), data = agg)
vc <- as.data.frame(VarCorr(m))
icc <- vc$vcov[vc$grp == "dataset"] / sum(vc$vcov)
cat(sprintf("\nICC(dataset) = %.3f  → %s\n", icc,
  ifelse(icc > 0.5, "배치 지배적: donor 분석에 dataset 랜덤효과 + within-dataset 메타분석 필수",
         "배치 영향 제한적: donor 수준 분석 진행 가능")))
print(summary(m)$coefficients)
print(confint(m, method = "Wald"))

# dataset 별 임상 구성 (배치가 생물학과 교란되는지 확인)
for (v in intersect(c("idh","grade","cell_type","tumor_type"), colnames(obj@meta.data)))
  print(table(as.character(obj[[ds_col]][,1]), as.character(obj[[v]][,1])))

write_result(agg, "donor_pseudobulk_SLC7A7.csv")
write_result(data.frame(ICC_dataset = icc, n_donor = nrow(agg),
                        n_dataset = dplyr::n_distinct(agg$dataset),
                        depth_beta = summary(m)$coefficients["scale(depth)", 1]),
             "ICC_gate_summary.csv")
message("완료: ICC 관문")
