source(here::here("R", "setup.R"))
suppressPackageStartupMessages({ library(Seurat); library(UCell); library(ggpubr) })

seu <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))

## 세포유형별 — donor 수준 요약
df <- data.frame(expr = FetchData(seu, vars = GENE_OF_INTEREST)[, 1],
                 ct = as.character(seu$cell_type), donor = seu$donor_id)


## AUCell (depth 강건) + donor 연속 검정
# 순환논리를 해결하기 위해 gene set에서 SLC7A7을 제거하고 분석진행
tr_set <- setdiff(build_transporter_geneset(universe = rownames(seu)), GENE_OF_INTEREST)
auc <- score_auc(seu, list(transport = tr_set,
                           arg_enzyme = build_arg_enzyme_geneset(universe = rownames(seu))))
seu$AUC_transport <- auc[, "transport"]; seu$AUC_argenzyme <- auc[, "arg_enzyme"]

r <- donor_continuous_test(seu, "AUC_transport")
print(r$coef); if (!is.null(r$mixed_coef)) print(r$mixed_coef)
cat(sprintf("within-dataset meta: beta=%.4g, se=%.4g, p=%.3g (%s datasets)\n",
            r$meta$pooled, r$meta$se, r$meta$p, r$meta$n_dataset))


# Estimate   Std. Error   t value     Pr(>|t|)
# (Intercept)  1.860909e-02 6.550550e-03  2.840844 0.0058085062
# gene        -4.998809e-03 3.871997e-03 -1.291016 0.2007168106
# depth        1.554954e-05 4.393760e-06  3.539006 0.0006982173
# Estimate   Std. Error       df    t value    Pr(>|t|)
# (Intercept) 1.836325e-02 6.607304e-03 73.61618 2.77923424 0.006911211
# gene        3.212899e-04 3.755128e-03 72.92202 0.08556032 0.932050427
# depth       1.386373e-05 4.410025e-06 73.08652 3.14368458 0.002411212



write_result(r$data, "Fig5_donor_AUC_transport_microglia.csv")
write_result(r$meta$per_dataset, "Fig5_within_dataset_meta_microglia.csv")

saveRDS(seu, file.path(DIR_DATA_PROC, "microglia_scored.rds"))
saveRDS(auc, file.path(DIR_DATA_PROC, "microglia_auc.rds"))



message("완료: Fig5")




"=================="

# add ---------------------------------------------------------------------

#12개 중 7개 양수, 5개 음수(Mathewson2021은 −0.614). 이건 "일관되게 재현됨"이 아니라 이질적입니다.
#최소 donor 수를 8~10으로 올리고 랜덤효과 메타 + I²로 다시 계산하세요.

per <- r$meta$per_dataset
per <- per[per$n >= 8, ]                       # 소표본 dataset 제외
w  <- 1/per$se^2; b <- sum(w*per$beta)/sum(w)
Q  <- sum(w*(per$beta - b)^2); df <- nrow(per) - 1
tau2 <- max(0, (Q - df)/(sum(w) - sum(w^2)/sum(w)))   # DerSimonian–Laird
wr <- 1/(per$se^2 + tau2); br <- sum(wr*per$beta)/sum(wr); ser <- sqrt(1/sum(wr))
c(beta = br, se = ser, p = 2*pnorm(-abs(br/ser)), I2 = max(0,(Q-df)/Q)*100, k = nrow(per))

#8 결과 
# beta           se            p           I2            k 
# 0.003226955  0.007511548  0.667487214 74.936371591  4.000000000
## 상당히 안 좋은 결과. study 간의 이질성이 74%다. 모든 그룹 간에서 공통된 효과가 아님.


#0 결과
# beta           se            p           I2            k 
# -0.002802993  0.006575844  0.669921318 62.891124756  9.000000000

