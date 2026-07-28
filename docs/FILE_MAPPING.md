# 원본 21파일 ↔ 재구성 구조 대응 / Original → refactored mapping

원본은 탐색·실행이 뒤섞인 스크립트였고, 재구성본은 로직을 **함수(R/)** 와
**실행 스크립트(analysis/)** 로 분리했습니다. 한 원본이 여러 곳으로 흩어지기도 합니다.

| 원본 파일 | 재구성 위치 |
|-----------|-------------|
| `TCGA.R` | `analysis/01`, `03` + `R/tcga_data.R`, `R/survival.R` |
| `TCGA_GBM.R` | `analysis/01`, `02` + `R/tcga_data.R`, `R/survival.R` |
| `egdeR.R` | `analysis/01` + `R/tcga_data.R`(normalize_logcpm, read_rawcounts_txt) |
| `TCGA_allgenes.R` | `analysis/02` + `R/survival.R`(safe_cox, analyze_gene_cox, plot_forest) |
| `TCGA_single.R` | `analysis/03` + `R/survival.R`(add_gene_expression, analyze_gene_single, plot_km) |
| `TCGA GSEA.R` | `analysis/04` + `R/gsea.R`(correlation_ranked_list, run_gsea_hallmark, run_gsego) |
| `GSEA_LGG.R` | `analysis/04` + `R/gsea.R`(run_gsego) |
| `GSEA_pseudobulk.R` | `analysis/08` + `R/gsea.R`(run_fgsea_pseudobulk, plot_nes_bar) |
| `GBmap 분석.R` | `analysis/05`, `06` + `R/scrna.R`(load_gbmap_seurat, find_slc_deg) |
| `patient_pseudo_bulk.R` | `analysis/07` + `R/scrna.R`(run_pseudobulk_limma) |
| `cellchat.R` | `analysis/09` + `R/scrna.R`(build_cellchat_cohort, run_cellchat, make_pos_neg_objects) |
| `CCC.R` | 초기 CellChat 탐색본 — 09로 대체(통합). 참고: 09가 최종 rigorous 버전 |
| `pseudotime.R` | `analysis/10` |
| `slingshot.R` | `analysis/11` |
| `volcano.R` | `analysis/12`(volcano/violin) |
| `seurat그래프그리는코드.R` | `R/utils.R`(theme_blank_frame) + `analysis/12` |
| `MR위한 SLC mapping 임시 R 파일.R` | `analysis/exploratory/01_mr_mapping.R` |
| `MR위한 cox_filtered mapping 임시 R 파일.R` | `analysis/exploratory/01_mr_mapping.R` (통합) |
| `MR.R` | `R/mr.R` + `exploratory/02` (초기 버전 → 통합본으로 대체) |
| `MR자동화.R` | `R/mr.R` + `exploratory/02` (중간 버전 → 통합) |
| `MRcML자동화.R` | `R/mr.R`(run_cml_for_gene) + `exploratory/02` |
| `MR통합 찐 최종 코드.R` | `R/mr.R` 전체 + `exploratory/02` (이 최종본 기준으로 정리) |

**중복 통합 메모**: TCGA-LGG/GBM 다운로드는 `download_tcga_expression(project)` 하나로,
LGG/GBM Cox는 cohort spec 인자로, MR 4개 버전은 최종 통합본 함수 세트로 합쳤습니다.
