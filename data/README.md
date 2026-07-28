# data/  (git 미추적 / not tracked)

## raw/ — 원본·외부 데이터 (직접 배치)
| 파일 | 설명 |
|------|------|
| `GBmap_core.h5ad` (바탕화면) | GBmap core 단일세포 atlas (h5ad) |

※ `{LGG,GBM}_TCGA_norm.csv`(logCPM)는 이제 `01_tcga_download_prepare.R` 이
  `normalize_logcpm()`(edgeR cpm)으로 **자동 생성**하여 processed/ 에 저장합니다.
  기존에 만들어 둔 norm csv 가 있으면 data/raw/ 에 넣어두면 그걸 우선 사용합니다.

## processed/ — 파이프라인 산출 중간물 (자동 생성)
| 파일 | 생성 | 사용처 |
|------|------|--------|
| `{LGG,GBM}_counts_clinical.rds` | 01 | 02, 03 |
| `{LGG,GBM}_download_info.txt` (다운로드 일시·샘플수 등) | 01 | 기록용 |
| `{LGG,GBM}_TCGA_rawcounts.txt` | 01 | (백업) |
| `{LGG,GBM}_TCGA_norm.csv` (logCPM) | 01 (edgeR cpm) | 02, 03, 04 |
| `{LGG,GBM}_all_cox_zscore.csv`  | 02 | forest, MR mapping |
| `seurat_gbmap.rds`, `microglia.rds` | 05 | 06–11 |
| `pb_fit.rds` (fit + paired_donors) | 07 | 08, 09 |
| `SLC_coxtest.csv` | exploratory/01 | exploratory/02 |
| `cds_microglia.rds` | 10 | (재사용) |

원본 코드가 참조하던 세션 이미지(.RData: seurat_0711, cellchat_0717, pseudotime_0720,
cds_0514 등)는 재구성본에서 명시적 .rds 산출물로 대체했습니다.
