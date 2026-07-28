# 파이프라인 · 데이터 흐름 / Pipeline data flow

## 실행 순서 (권장)
```
[bulk 먼저 — depth 문제 없음, 2~3주면 논문 절반]
01 TCGA 준비 ─▶ 02 gene set ─▶ 03 LASSO/Cox(Fig2) ─▶ 04 다중코호트 검증(Fig3) ─▶ 05 GSEA/decon(Fig4)

[그 다음 단일세포 — 00 관문 통과 시]
06 GBmap 전처리 ─▶ 00 ICC 관문 ─▶ 07 atlas(Fig5) ─▶ 08 myeloid(Fig6) ─▶ 09 trajectory/CCC(Fig7)

[마지막]
10 공간(Fig8) · 11 depth 민감도(Supp) · 12 그림 틀
```

## A. Bulk 축
```
01_tcga_download_prepare.R
   get_tcga_data(cache) → prepare_count_matrix → normalize_logcpm(edgeR)
     → data/processed/{LGG,GBM}_counts_clinical.rds, {LGG,GBM}_TCGA_norm.csv

02_curate_transporter_geneset.R
   build_transporter_geneset()  = MSigDB 3종 ∪ 수동 y+L/CAT
     → results/transporter_geneset.csv        ★ Fig2 입력

03_lasso_cox_signature.R   [Fig2]
   run_lasso_cox(10-fold CV) → selected genes → compute_risk_score()
   univariate_cox() / multivariable_cox() / stratified_cox(by=idh)
     → results/LASSO_selected_genes.csv, multivariable_cox_train.csv
   ※ SLC7A7 탈락 시 candidate-gene 설계로 전환(코드가 양쪽 산출)

04_multicohort_validation.R [Fig3]
   load_all_cohorts()  (TCGA_LGG/GBM, CGGA_693, CGGA_325 …)
   km_cohort() · timedep_roc(12/36/60mo) · dca_cohort() · multivariable_cox()
     → results/Fig3_*.csv, figures/Fig3_*.pdf

05_bulk_clinical_gsea_decon.R [Fig4]
   임상 violin · correlation_ranked_list→GSEA · run_estimate()/run_xcell()
   correlate_gene_fractions() · adjust_for_fraction()   ← myeloid 보정 후 예후 유지 여부
```

## B. 단일세포 축 (donor 수준)
```
06_sc_preprocess_gbmap.R
   load_gbmap_bpcells()  BPCells 온디스크 + obs 전체 + UMAP + LogNormalize
     → seurat_gbmap.rds, microglia.rds, mac.rds, myeloid.rds

00_donor_pseudobulk_icc.R   ★관문
   lmer(expr ~ depth + (1|dataset)) → ICC(dataset)
   ICC 높으면 이후 donor 분석에서 dataset 을 반드시 보정

07_sc_atlas_donor.R [Fig5]  score_auc(AUCell) → donor_continuous_test(+depth)
                            SLC3A2 공발현 donor 회귀
08_sc_myeloid_donor.R [Fig6] donor_subtype_proportion() ← 핵심
                            compare_mg_mdm(), run_pseudobulk_limma → fgsea
                             → data/processed/pb_fit.rds
09_trajectory_cellchat.R [Fig7] monocle3(연속값 투영) · slingshot · CellChat(donor high/low)
```

## C. 공간 / 민감도
```
10_spatial.R [Fig8]  load_spatial_all → add_spatial_modules(hypoxia/myeloid/necrosis/transport)
                     spatial_spot_correlation(SLC7A7 vs ASS1/ASL/SLC3A2) · spatial_niche_test
11_depth_sensitivity.R [Supp] depth_sensitivity(): depth비 · glmer OR · 매칭 음성대조
```

## 공유 객체 / 모듈 함수
| 모듈 | 주요 함수 |
|------|-----------|
| `R/gene_sets.R` | `build_transporter_geneset`, `build_arg_enzyme_geneset`, `get_msig_genes` |
| `R/cohorts.R` | `load_cohort`, `load_all_cohorts`, `harmonise_clinical`, `cohort_gene` |
| `R/signature.R` | `run_lasso_cox`, `compute_risk_score`, `univariate_cox`, `multivariable_cox`, `stratified_cox` |
| `R/validation.R` | `km_cohort`, `timedep_roc`, `dca_cohort`, `forest_multicohort` |
| `R/deconvolution.R` | `run_estimate`, `run_xcell`, `correlate_gene_fractions`, `adjust_for_fraction` |
| `R/sc_donor.R` | `donor_pseudobulk`, `score_auc`, `donor_continuous_test`, `donor_subtype_proportion`, `depth_sensitivity`, `compare_mg_mdm` |
| `R/spatial.R` | `load_spatial_all`, `add_spatial_modules`, `spatial_spot_correlation`, `spatial_niche_test` |
| `R/scrna.R` | `load_gbmap_bpcells`, `read_obs_all`, `run_pseudobulk_limma`, `run_cellchat` |
