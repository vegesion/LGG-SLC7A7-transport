# SLC7A7 and microglia in Low-Grade Glioma / Glioblastoma

> **한국어 + English** — 뇌종양(LGG·GBM)에서 **SLC7A7 발현 양성/음성(+/–) microglia**가
> 종양 미세환경에 미치는 영향을 bulk RNA-seq, 단일세포(scRNA-seq), Mendelian
> Randomization으로 분석한 재현 가능한 파이프라인입니다.

이 저장소는 흩어져 있던 21개 개별 R 스크립트를, 실제 연구 랩에서 쓰는 방식대로
**설정(config) · 재사용 함수(R/) · 번호형 실행 파이프라인(analysis/)** 구조로 재구성한 것입니다.
분석 절차와 통계(표준화·임계값·시드·검정법)는 원본과 동일하게 보존했습니다.

---

## 핵심 설계 / Design

핵심은 **"변수·함수의 연동"** 입니다. 예전에는 파일마다 변수를 그때그때 바꿔가며
같은 코드를 복붙했지만, 이제는:

- **`config/config.R`** — 모든 경로·파라미터·임계값·대상 유전자(`GENE_OF_INTEREST = "SLC7A7"`)를 한 곳에서 관리. 절대경로 하드코딩 없음(`here::here()` 상대경로).
- **`R/`** — 재사용 함수 라이브러리. 예: `safe_cox()`, `run_cox_screen()`, `run_pseudobulk_limma()`, `run_cellchat()`, `run_mr_for_gene()`. LGG/GBM처럼 반복되던 코드는 인자(cohort spec)만 바꿔 하나의 함수로 처리.
- **`analysis/`** — 얇은 실행 스크립트. 각 스크립트는 `source(here::here("R","setup.R"))` 한 줄로 전체 환경을 불러온 뒤, 함수를 호출해 결과를 `results/`·`figures/`에 저장.

모든 스크립트는 첫 줄에 **목적·입력·출력·원본파일**을 주석으로 명시합니다.

---

## 실행 방법 / How to run

1. RStudio에서 **`LGG-SLC7A7-microglia.Rproj`** 를 엽니다 (작업 디렉터리가 루트로 고정).
2. 최초 1회 의존성 설치: `source("install_dependencies.R")`
3. OpenGWAS 토큰 설정(МR용): `.Renviron.example` → `.Renviron` 복사 후 토큰 입력.
4. 원본 데이터를 `data/raw/`에 배치 (아래 데이터 표 참고).
5. `analysis/` 스크립트를 번호 순서(01→12, 그다음 exploratory)로 실행.

```r
source(here::here("R", "setup.R"))     # 모든 스크립트 공통 첫 줄
```

---

## 파이프라인 개요 / Pipeline

**전략**: 선행연구(ASL, 아르기닌 *대사효소* 축)와 같은 틀·다른 축.
본 연구는 **수송 축(y+L/CAT)** · **donor 수준 연속 분석** · **공간전사체** · **IDH 층화 후 독립성**으로 차별화한다.
**Fig 2–4(bulk)를 먼저 완성**하면 그 자체로 논문이 성립하고, 단일세포·공간은 그 위에 얹는 층이다.

| # | 스크립트 | Figure | 내용 |
|---|----------|--------|------|
| 00 | `00_donor_pseudobulk_icc.R` | — (관문) | dataset ICC: 배치 vs donor 분산 분해 |
| 01 | `01_tcga_download_prepare.R` | — | TCGA 다운로드 + logCPM(edgeR) |
| 02 | `02_curate_transporter_geneset.R` | Fig 2 준비 | 수송체 gene set curation (MSigDB 3종 + y+L/CAT 수동) |
| 03 | `03_lasso_cox_signature.R` | **Fig 2** | LASSO Cox + 10-fold CV → 유전자 선정, 단/다변량 Cox |
| 04 | `04_multicohort_validation.R` | **Fig 3** | TCGA/CGGA 693+325 KM · time-dep ROC · DCA · **IDH 보정/층화** |
| 05 | `05_bulk_clinical_gsea_decon.R` | **Fig 4** | 임상 상관 · GSEA · **면역 deconvolution + myeloid 보정** |
| 06 | `06_sc_preprocess_gbmap.R` | — | GBmap → Seurat v5 (BPCells 온디스크) |
| 07 | `07_sc_atlas_donor.R` | **Fig 5** | atlas UMAP · AUCell(수송) · SLC3A2 공발현(depth 보정) |
| 08 | `08_sc_myeloid_donor.R` | **Fig 6** | **donor 수준 발현 vs 아종 비율 회귀** · MG vs MDM · pseudobulk DEG/GSEA |
| 09 | `09_trajectory_cellchat.R` | **Fig 7** | Monocle3(**연속값 투영**) · Slingshot · CellChat(donor high/low) |
| 10 | `10_spatial.R` | **Fig 8** | 공간 분포 · niche · **spot 수준 SLC7A7 vs ASS1** |
| 11 | `11_depth_sensitivity.R` | Supp | depth 비 · glmer OR · **매칭 음성대조** |
| 12 | `12_figures.R` | — | 논문용 "그림 틀"(PPT 편집) |
| S | `supplementary/celllevel_slc7a7_deg.R` | Supp | 세포 단위 +/- (main 사용 금지) |
| E | `exploratory/MR/*` | — | MR (논문 미사용) |

### 설계 원칙 (이 프로젝트의 방법론적 차별점)
1. **분석 단위는 donor** — 33만 세포 단위 검정은 pseudoreplication 으로 p 가 무의미해진다.
2. **발현은 연속값** — 세포 단위 +/- 이분은 sequencing depth 와 교란되므로 main figure 에서 제외.
3. **depth 보정 내장** — donor 수준 회귀에 `nFeature_RNA` 를 공변량으로, gene set 점수는
   순위 기반 **AUCell**(depth 강건) 사용.
4. **모듈 점수에서 self-gene 제외** — SLC7A7 로 그룹을 나눈 뒤 SLC7A7 이 포함된 모듈을
   비교하는 순환논리를 차단.

스크립트 간 데이터 흐름(공유 객체·중간 CSV)은 [`docs/PIPELINE.md`](docs/PIPELINE.md) 참고.
원본 21파일 ↔ 새 구조 대응은 [`docs/FILE_MAPPING.md`](docs/FILE_MAPPING.md).

---

## "그림 틀만" 워크플로 / Figure-frame workflow

`12_figures.R`와 `theme_blank_frame()`은 축·제목·모든 텍스트를 제거하고 tick만 남긴
그림을 출력합니다. 이 벡터 그림(EMF/TIFF)을 PPT로 가져가 라벨을 다시 그려
논문 제출본을 만드는 기존 워크플로를 그대로 함수화한 것입니다.

---

## ⚠️ 보안 / Security

원본 MR 코드에는 OpenGWAS JWT 토큰이 하드코딩되어 있었습니다. 재구성본은 토큰을
**코드에서 제거**하고 `Sys.getenv("OPENGWAS_JWT")`(→ `.Renviron`)로만 읽습니다.
`.Renviron`·`data/`·`*.RData`는 `.gitignore`로 커밋에서 제외됩니다.
**공개 저장소로 올리기 전, 과거 커밋에 토큰이 남아있지 않은지 확인하세요.**

## 주의 / Caveat

이 저장소는 원본 로직을 보존한 구조 재구성입니다. 데이터·R 환경 없이 실행 검증은
불가능하므로, 처음 돌릴 때는 스크립트를 번호 순으로 한 단계씩 실행하며 확인하세요.

---

## 재현환경 / Reproducible environment (renv)

패키지 버전 고정은 **renv 로 본인이 직접 초기화**합니다. (임의로 만든 renv.lock 을 넣지
않았습니다 — lock 은 본인 컴퓨터의 실제 설치 상태에서 생성돼야 정확하기 때문입니다.)

```r
source("install_dependencies.R")   # (권장) 먼저 패키지 설치
install.packages("renv")
renv::init()                       # 코드 스캔 → 의존성 파악 → renv.lock 생성
```

이후 `renv::snapshot()`(lock 갱신) / `renv::restore()`(복원)로 관리합니다.
자세한 절차·문제 해결은 [`docs/RENV_SETUP.md`](docs/RENV_SETUP.md) 를 보세요.
renv 는 **선택 사항**이며, 패키지만 설치돼 있으면 renv 없이도 파이프라인은 실행됩니다.

---

## 데이터 경로 설정 / Data paths

큰 원본 파일(h5ad, 정규화 csv)은 저장소로 복사하지 않고 **기존 폴더를 직접 참조**합니다.
`config/config.R` 의 `DATA_SEARCH_ROOTS` 가 후보 폴더들을 순서대로 탐색합니다:

1. `data/raw/` (저장소 내부, 있으면 최우선)
2. `C:/Users/helis/바탕화면` — **GBmap h5ad 위치(확인됨)**: `51c6f87a-…-732b4.h5ad`
3. OneDrive `…/bioinformatics/analysis/files`, `…/analysis`, `…/analysis/LGG-SLC4A7`
4. G드라이브 백업 경로

정규화 csv(`LGG_TCGA_norm.csv`, `GBM_TCGA_norm.csv`) 위치는 아직 미확인 상태라,
원본 코드 기준 후보 경로를 넣어 두었습니다. 실제 위치를 알게 되면
`config/config.R` 의 `DATA_SEARCH_ROOTS` **맨 위에 한 줄만 추가**하면 됩니다:

```r
DATA_SEARCH_ROOTS <- c(
  "D:/내가_찾은_데이터_폴더",   # ← 여기에 실제 경로 추가
  DIR_DATA_RAW,
  ...
)
```

파일을 못 찾으면 스크립트가 `require_data()` 를 통해 **어떤 파일을 어디서 찾았는지**
명확히 알려주며 멈춥니다. 경로 구분자는 `/` 또는 `\\` (역슬래시 2개) 를 쓰세요.

---

## 문제 해결 / Troubleshooting

**`Error in readRDS(dest): 알 수 없는 입력 포맷 (unknown input format)`**
저장소 메타데이터(PACKAGES.rds)나 패키지 파일 다운로드가 깨졌다는 뜻입니다
(받다 만 파일, 또는 프록시/VPN이 HTML 에러 페이지를 반환한 경우). 해결:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"),
        download.file.method = "libcurl", timeout = 600)
install.packages("BiocManager")
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
```

그 뒤 `source("install_dependencies.R")` 재실행. 이 설정은 스크립트에도 반영돼 있습니다.
그래도 나면: ① **새 R 세션**에서 재시도, ② 사내 프록시/VPN 잠시 해제, ③ 다른 CRAN 미러 사용.

**`... replaces Bioconductor standard repositories` 경고** — 무해합니다.
`BiocManager::install(version = "3.20")` 을 한 번 실행하면 저장소가 올바로 등록됩니다.

**`node stack overflow` (renv 실행 중)**
이전에 남은 `renv/` 폴더나 `.Rprofile` 과 충돌해 `renv::init()` 의 세션 재시작이
무한 재귀에 빠지는 경우입니다. 완전 초기화 후 다시 시도하세요:

```r
unlink("renv", recursive = TRUE)
if (file.exists(".Rprofile")) file.rename(".Rprofile", ".Rprofile.bak")
# → R 재시작(Ctrl+Shift+F10) 후
renv::init(restart = FALSE)
```

세션이 아예 안 열리면 탐색기에서 `.Rprofile` 이름을 바꾼 뒤 R 을 재시작하세요.
자세한 내용은 [`docs/RENV_SETUP.md`](docs/RENV_SETUP.md).

**TCGA 다운로드가 반복 실패 / 중간에 잘림(`chunks ... was not correct`)**
청크가 실제 크기보다 훨씬 작은 지점(예: 11MB)에서 "completed" 되며 반복 실패하면,
백신·프록시가 다운로드를 끊는 **truncation** 입니다. 순서대로:

```r
# 1) 청크를 최소로 (config/config.R): TCGA_FILES_PER_CHUNK <- 1
options(timeout = 10000, download.file.method = "libcurl")
unlink(DIR_GDC, recursive = TRUE)
source(here::here("analysis", "01_tcga_download_prepare.R"))
```
+ 백신 실시간검사 잠시 끄기 / 다른 네트워크(핫스팟). 그래도 안 되면 **gdc-client**
(`TCGA_DOWNLOAD_METHOD <- "client"`)로 전환하세요. 자세한 절차: [`docs/GDC_DOWNLOAD.md`](docs/GDC_DOWNLOAD.md).

**h5ad 로드 시 RAM 폭발 (32GB RAM 인데도 7GB 파일에서 터짐)**
`anndata::read_h5ad()` (R 래퍼)는 `backed = "r"` 로 열어도 `ad$X` 에 접근하는 순간
**전체를 dense R 행렬로 변환**해 메모리를 폭발시킵니다. → **anndata 래퍼를 쓰지 마세요.**

대신 (1) 크기를 먼저 재고 → (2) hdf5r 로 sparse 그대로 읽습니다(파이썬 미경유):

```r
source(here::here("R", "setup.R"))
h5ad_info(PATH_H5AD)                 # 데이터 미읽음: cells/genes/nnz + 예상 메모리(GB)
dat <- read_h5ad_sparse(PATH_H5AD)   # sparse dgCMatrix (dense 변환 없음)
print(format(object.size(dat$counts), units = "GB"))
```

`h5ad_info()` 의 "예상 메모리" 가 대략 20GB 이하면 32GB 에서 안전하게 읽힙니다.
`load_gbmap_seurat()` 도 기본이 이 hdf5r 경로예요. (예전 방식은 `use_python = TRUE`.)

**그래도 너무 커서 안 들어가면 (예상 메모리 > ~20GB)** — 온디스크 방식 **BPCells** 를 쓰세요.
데이터를 RAM 에 다 올리지 않고 디스크에 둔 채 Seurat v5 로 분석합니다(대형 atlas 표준 해법):

```r
remotes::install_github("bnprks/BPCells/r")
library(BPCells)
# h5ad 의 X 를 온디스크 행렬로 1회 변환 후 Seurat 로 사용 (자세히는 BPCells 문서)
```

원인은 reticulate/anndata 버전이 sparse→dense 로 바뀐 것이며, hdf5r 직접 읽기는
그 버전 의존성 자체를 없앱니다.

**`source("R/setup.R")` 한 줄인데 몇 분씩 멈춤 (hang)**
config 가 데이터를 찾으려 후보 경로에 `file.exists()` 를 도는데, 그 목록에 **클라우드/네트워크
경로**(OneDrive 온라인전용, Google Drive `G:/다른 컴퓨터`)가 있으면 각 확인이 네트워크 응답을
기다리며 멈춥니다. → `config/config.R` 의 `DATA_SEARCH_ROOTS` 를 **로컬 경로만** 남기세요
(기본값이 이미 로컬 전용으로 바뀌었습니다). 멈춘 R 은 컴퓨터를 끄지 말고 **Esc / Session→Restart R
(Ctrl+Shift+F10)** 로만 중단하세요.

**대용량 작업 시 RStudio 설정 (강제종료·워크스페이스 손상 예방)**
7GB h5ad 같은 큰 객체를 다루면, RStudio 가 종료 때 `.RData` 로 저장하고 시작 때 복원하려다
멈추거나 세션이 깨질 수 있습니다. **Tools → Global Options → General** 에서:
- *Save workspace to .RData on exit* → **Never**
- *Restore .RData into workspace at startup* → **체크 해제**
로 두세요. (분석 결과는 스크립트가 명시적으로 `saveRDS()` 로 저장하므로 워크스페이스 자동저장은 불필요)
