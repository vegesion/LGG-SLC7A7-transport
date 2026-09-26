# =============================================================================
# 08c_sc_donor_gsea_heatmap.R      [Fig 6 보조 — donor별 GSEA 일관성 히트맵]
# -----------------------------------------------------------------------------
# 목적 (역할 분담을 반드시 지킬 것)
#   · 통계적 주장  = 08b 의 within/between 혼합모형 GSEA  (본문 표 / main figure)
#   · 이 스크립트  = "그 결과가 소수 donor 의 인공물이 아니라 다수 donor 에
#                    걸친 실체다" 를 눈으로 확인시키는 **시각화**
#   → 히트맵의 셀은 donor 내부 '세포 수준' GSEA 이므로 pseudoreplication 을
#     피하지 못한다. 이 그림에서 나온 p 값은 논문의 근거로 쓰지 않는다.
#     (섹션 5 의 sign test 도 vote-counting 보조 지표일 뿐이다.)
#
# 설계
#   donor d 안에서 유전자별 OLS:  gene ~ SLC7A7 + depth
#     → t(β_SLC7A7) 로 랭킹 → fgsea(Hallmark 50)  → NES_{pathway, donor}
#   이는 08b 혼합모형의 g_within 항과 정확히 같은 질문을 donor 별로 분해한 것
#   (donor 내부만 쓰므로 donor 중심화가 자동으로 이루어짐).
#
#   ★ depth 음성대조: 같은 적합에서 나온 t(β_depth) 로 동일 GSEA 를 한 번 더
#     돌린다. IFN 계열이 depth 랭킹에서도 똑같이 켜지면 이 그림 전체가
#     depth 교란이라는 뜻이다. 반드시 함께 보고할 것.
#
# 입력 : data/processed/microglia.rds
# 출력 : results/Fig6_donor_gsea_long.csv          (donor × pathway 전체)
#        results/Fig6_donor_gsea_qc.csv            (donor 별 포함/제외 사유)
#        results/Fig6_donor_gsea_consistency.csv   (pathway 별 일관성 요약)
#        figures/Fig6_donor_gsea_heatmap.pdf       (본 그림)
#        figures/Fig6_donor_gsea_heatmap_depth.pdf (음성대조)
#        figures/Fig6_donor_gsea_consistency.pdf
#        data/processed/donor_gsea_raw.rds         (재렌더링용 중간 산물)
# =============================================================================

source(here::here("R", "setup.R"))
suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
  library(fgsea); library(msigdbr)
  library(dplyr); library(tidyr); library(ggplot2)
})

## ── 0. 설정 ─────────────────────────────────────────────────────────────────
MIN_CELLS_DONOR   <- 100    # ① donor 당 최소 세포 수  ← 요청값
MIN_POS_CELLS     <- 20     # ② donor 당 최소 SLC7A7+ 세포 수   ★ 아래 주석 참조
MIN_DETECT_RATE   <- 0.03   # ③ donor 내 SLC7A7 검출률 하한
GENE_MIN_PCT      <- 0.10   # 전역 유전자 필터 (풀링 기준 검출률)
GENE_MIN_CELLS_D  <- 3      # donor 내 최소 검출 세포 수 (0분산 유전자 제거)
NPERM_SIMPLE      <- 10000  # donor 수만큼 반복하므로 08b(1e5)보다 낮춤
SEED              <- 42

ROW_ORDER   <- "median_nes"   # "median_nes" | "frac_same"
COL_ORDER   <- "anchor"       # "anchor" | "concordance" | "ncells" | "slc7a7"
ANCHOR_PATH <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
GREY_NONSIG <- TRUE           # donor 내 BH FDR ≥ 0.05 셀을 회색 처리
SIG_CUT     <- 0.05
FACET_BY_DATASET <- FALSE     # TRUE 로 두면 배치 교란 점검용 버전이 나옴
RENDER_ONLY <- FALSE          # TRUE 면 저장된 RDS 로 그림만 다시 그림

# ★ MIN_POS_CELLS / MIN_DETECT_RATE 를 넣은 이유 (이게 이 스크립트의 핵심 방어선)
#   세포 수 100 을 넘겨도 그 donor 에서 SLC7A7 이 3개 세포에만 잡히면,
#   β_SLC7A7 는 사실상 그 3개 세포가 결정한다. NES 는 난수가 되지만
#   "세포 100개 이상"이라는 기준은 통과해 버린다.
#   즉 GSEA 안정성을 결정하는 건 총 세포 수가 아니라 **SLC7A7 정보량**이다.
#   두 관문의 탈락 donor 수를 각각 출력하니, 어느 쪽이 실제 제약인지 확인할 것.

set.seed(SEED)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

RAW_PATH <- file.path(DIR_DATA_PROC, "donor_gsea_raw.rds")

# =============================================================================
# 1. 데이터 로드 · donor 관문
# =============================================================================
if (!RENDER_ONLY) {

mg <- readRDS(file.path(DIR_DATA_PROC, "microglia.rds"))

# ⚠ 객체 이름이 microglia.rds 라도 실제 내용이 TAM 전체(MG + BDM)라면
#   figure legend·본문에서 "microglia" 로 쓰면 안 된다(금지 프레이밍).
#   아래 출력으로 실제 구성을 확인하고 캡션 문구를 결정할 것.
ct_col <- intersect(c("annotation_level_3", "annotation_level_2", "cell_type"),
                    colnames(mg@meta.data))[1]
if (!is.na(ct_col)) {
  cat("[구성] ", ct_col, " 분포\n", sep = "")
  print(sort(table(mg@meta.data[[ct_col]]), decreasing = TRUE))
}

batch_col <- find_batch_col(mg)
if (is.na(batch_col)) stop("dataset(배치) 컬럼을 찾지 못했습니다.")

meta <- data.frame(
  cell    = colnames(mg),
  g       = Seurat::FetchData(mg, vars = GENE_OF_INTEREST)[, 1],
  depth   = mg[[DEPTH_COVARIATE]][, 1],
  donor   = as.character(mg$donor_id),
  dataset = as.character(mg[[batch_col]][, 1]),
  stringsAsFactors = FALSE)

## donor 별 QC 지표
qc <- meta %>%
  group_by(donor, dataset) %>%
  summarise(
    n_cells    = dplyr::n(),
    n_pos      = sum(g > 0),
    detect     = mean(g > 0),
    sd_g       = stats::sd(g),
    mean_g     = mean(g),
    cor_g_depth = suppressWarnings(stats::cor(g, depth, method = "spearman")),
    .groups = "drop") %>%
  mutate(
    fail_cells  = n_cells < MIN_CELLS_DONOR,
    fail_pos    = n_pos   < MIN_POS_CELLS,
    fail_detect = detect  < MIN_DETECT_RATE,
    fail_var    = !is.finite(sd_g) | sd_g == 0,
    keep = !(fail_cells | fail_pos | fail_detect | fail_var))

cat(sprintf(
"\n[관문] 전체 donor %d\n  ① 세포수 < %d          탈락 %d\n  ② SLC7A7+ 세포 < %d    탈락 %d\n  ③ 검출률 < %.0f%%        탈락 %d\n  ④ 분산 0               탈락 %d\n  → 최종 포함 donor %d\n",
  nrow(qc), MIN_CELLS_DONOR, sum(qc$fail_cells),
  MIN_POS_CELLS, sum(qc$fail_pos),
  MIN_DETECT_RATE * 100, sum(qc$fail_detect),
  sum(qc$fail_var), sum(qc$keep)))

if (any(abs(qc$cor_g_depth[qc$keep]) > 0.5, na.rm = TRUE))
  warning(sprintf("SLC7A7–depth |rho| > 0.5 인 donor %d명 — 해당 열은 해석 주의",
                  sum(abs(qc$cor_g_depth[qc$keep]) > 0.5, na.rm = TRUE)))

write_result(qc, "Fig6_donor_gsea_qc.csv")

keep_donors <- qc$donor[qc$keep]
if (length(keep_donors) < 10)
  stop("포함 donor 가 10명 미만입니다. 관문을 재검토하십시오.")
meta <- meta[meta$donor %in% keep_donors, ]

# =============================================================================
# 2. 전역 유전자 universe (donor 간 NES 비교 가능성을 위해 고정)
# =============================================================================
expr <- Seurat::GetAssayData(mg, assay = "RNA", layer = "data")
expr <- as(expr[, meta$cell, drop = FALSE], "dgCMatrix")   # genes × cells
pct   <- Matrix::rowMeans(expr > 0)
genes <- setdiff(rownames(expr)[pct >= GENE_MIN_PCT], GENE_OF_INTEREST)  # self 제외
expr  <- expr[genes, , drop = FALSE]
rm(mg); gc(verbose = FALSE)

message(sprintf("[유전자] universe %s개 (풀링 검출률 ≥ %.0f%%, self 제외)",
                format(length(genes), big.mark = ","), GENE_MIN_PCT * 100))

# =============================================================================
# 3. donor별 벡터화 OLS  →  t(β_SLC7A7), t(β_depth)
#    lmer/lm 반복 대신 정규방정식을 한 번에 푼다 (유전자 수천 개 × donor 수십 명).
# =============================================================================
donor_tstats <- function(dn) {
  idx   <- which(meta$donor == dn)
  n     <- length(idx)
  Yd    <- Matrix::t(expr[, idx, drop = FALSE])            # cells × genes
  # donor 내 0분산/저검출 유전자 제거
  det_n <- Matrix::colSums(Yd > 0)
  keepg <- det_n >= GENE_MIN_CELLS_D & det_n <= (n - 1)
  if (sum(keepg) < 500) return(NULL)
  Yd <- Yd[, keepg, drop = FALSE]

  X <- cbind(1,
             as.numeric(scale(meta$g[idx])),
             as.numeric(scale(meta$depth[idx])))
  colnames(X) <- c("int", "g", "depth")
  if (any(!is.finite(X))) return(NULL)
  XtXi <- try(solve(crossprod(X)), silent = TRUE)          # 공선성 방어
  if (inherits(XtXi, "try-error")) return(NULL)

  XtY <- as.matrix(Matrix::crossprod(X, Yd))               # 3 × G
  B   <- XtXi %*% XtY
  rss <- pmax(Matrix::colSums(Yd^2) - colSums(B * XtY), .Machine$double.eps)
  s2  <- rss / (n - ncol(X))

  t_g <- B["g", ]     / sqrt(s2 * XtXi["g", "g"])
  t_d <- B["depth", ] / sqrt(s2 * XtXi["depth", "depth"])
  nm  <- colnames(Yd)
  list(t_gene  = setNames(t_g[is.finite(t_g)], nm[is.finite(t_g)]),
       t_depth = setNames(t_d[is.finite(t_d)], nm[is.finite(t_d)]),
       n_cells = n, n_genes = sum(keepg))
}

HALLMARK <- get_hallmark_list()

run_fgsea_one <- function(stat_vec, dn, level) {
  sv <- sort(stat_vec[is.finite(stat_vec)], decreasing = TRUE)
  if (length(sv) < 500) return(NULL)
  set.seed(SEED)
  fg <- suppressWarnings(
    fgsea::fgsea(pathways = HALLMARK, stats = sv,
                 minSize = 10, maxSize = 500, eps = 0,
                 nPermSimple = NPERM_SIMPLE))
  if (is.null(fg) || nrow(fg) == 0) return(NULL)
  fg %>%
    dplyr::select(pathway, NES, pval, padj, size) %>%
    dplyr::mutate(donor = dn, level = level)
}

t0 <- Sys.time()
res <- pbapply::pblapply(keep_donors, function(dn) {
  ts <- donor_tstats(dn)
  if (is.null(ts)) return(NULL)
  list(gene  = run_fgsea_one(ts$t_gene,  dn, "SLC7A7"),
       depth = run_fgsea_one(ts$t_depth, dn, "depth"),
       info  = data.frame(donor = dn, n_cells = ts$n_cells, n_genes = ts$n_genes))
})
names(res) <- keep_donors
message(sprintf("[GSEA] donor %d명 · 소요 %.1f 분",
                sum(!vapply(res, is.null, logical(1))),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

gsea_long <- bind_rows(lapply(res, function(x) if (is.null(x)) NULL else
                              bind_rows(x$gene, x$depth)))
info_df   <- bind_rows(lapply(res, function(x) if (is.null(x)) NULL else x$info))

gsea_long <- gsea_long %>%
  left_join(qc %>% select(donor, dataset, n_pos, detect, mean_g, cor_g_depth),
            by = "donor") %>%
  left_join(info_df, by = "donor")

saveRDS(list(gsea_long = gsea_long, qc = qc), RAW_PATH)
write_result(gsea_long, "Fig6_donor_gsea_long.csv")

} else {
  obj <- readRDS(RAW_PATH); gsea_long <- obj$gsea_long; qc <- obj$qc
}

# =============================================================================
# 4. pathway 일관성 요약 (행 정렬 기준)
# =============================================================================
consist <- gsea_long %>%
  filter(level == "SLC7A7") %>%
  group_by(pathway) %>%
  summarise(
    n_donor     = dplyr::n(),
    median_nes  = median(NES, na.rm = TRUE),
    iqr_nes     = IQR(NES, na.rm = TRUE),
    n_pos       = sum(NES > 0, na.rm = TRUE),
    frac_same   = max(mean(NES > 0, na.rm = TRUE), mean(NES < 0, na.rm = TRUE)),
    dir_major   = ifelse(mean(NES > 0, na.rm = TRUE) >= 0.5, "up", "down"),
    n_sig       = sum(padj < SIG_CUT, na.rm = TRUE),
    frac_sig    = mean(padj < SIG_CUT, na.rm = TRUE),
    # ▼ vote-counting sign test: donor 1명 = 관측 1개 이므로 pseudoreplication 은
    #   아니지만, donor 별 NES 정밀도 차이(세포 100 vs 17,000)를 무시한다.
    #   ★ 보조 지표. 본문의 통계적 주장은 08b 혼합모형 GSEA 로만 한다.
    sign_p      = stats::binom.test(max(n_pos, n_donor - n_pos), n_donor, 0.5)$p.value,
    .groups = "drop") %>%
  mutate(sign_fdr = p.adjust(sign_p, "BH")) %>%
  arrange(desc(median_nes))
write_result(consist, "Fig6_donor_gsea_consistency.csv")

cat("\n===== 일관되게 양수 (상위 8) =====\n")
print(as.data.frame(head(consist, 8)), digits = 3)
cat("\n===== 일관되게 음수 (하위 8) =====\n")
print(as.data.frame(tail(consist, 8)), digits = 3)

# depth 음성대조와의 대조표 — 같은 pathway 가 양쪽에서 같은 방향이면 위험 신호
depth_consist <- gsea_long %>% filter(level == "depth") %>%
  group_by(pathway) %>%
  summarise(median_nes_depth = median(NES, na.rm = TRUE),
            frac_sig_depth   = mean(padj < SIG_CUT, na.rm = TRUE), .groups = "drop")
cmp <- consist %>% select(pathway, median_nes, frac_same, frac_sig) %>%
  left_join(depth_consist, by = "pathway") %>%
  mutate(same_sign_as_depth = sign(median_nes) == sign(median_nes_depth)) %>%
  arrange(desc(abs(median_nes)))
write_result(cmp, "Fig6_donor_gsea_vs_depth.csv")
cat(sprintf("\n[음성대조] SLC7A7 랭킹과 depth 랭킹의 median NES 부호가 같은 pathway: %d / %d\n",
            sum(cmp$same_sign_as_depth, na.rm = TRUE), nrow(cmp)))
cat(sprintf("           |median NES| 상위 10 중 부호 일치: %d\n",
            sum(head(cmp$same_sign_as_depth, 10), na.rm = TRUE)))

# =============================================================================
# 5. 히트맵
# =============================================================================
make_heatmap <- function(lv, ttl, sub) {
  
  d <- gsea_long %>% filter(level == lv)
  
  ## -------------------------------------------------------------------------
  ## 0. 모든 pathway가 non-significant인 donor 제거
  ##    - GSEA 결과 자체에서는 donor를 제거하지 않음
  ##    - heatmap 시각화에서만 제거
  ## -------------------------------------------------------------------------
  donor_keep <- d %>%
    group_by(donor) %>%
    summarise(
      has_sig = any(padj < SIG_CUT, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(has_sig) %>%
    pull(donor)
  
  d <- d %>%
    filter(donor %in% donor_keep)
  
  message(
    lv, ": ",
    length(donor_keep),
    " donors retained; ",
    length(setdiff(unique(gsea_long$donor[gsea_long$level == lv]),
                   donor_keep)),
    " all-nonsignificant donors removed."
  )
  
  ## 행 정렬 -----------------------------------------------------------------
  row_key <- if (ROW_ORDER == "frac_same") {
    consist %>%
      mutate(k = ifelse(dir_major == "up", frac_same, -frac_same)) %>%
      select(pathway, k)
  } else {
    consist %>%
      mutate(k = median_nes) %>%
      select(pathway, k)
  }
  
  row_lvl <- row_key %>%
    arrange(k) %>%
    pull(pathway)   # 아래→위 = 음수→양수
  
  
  ## 열 정렬 -----------------------------------------------------------------
  col_key <- switch(
    COL_ORDER,
    
    "anchor" = d %>%
      filter(pathway == ANCHOR_PATH) %>%
      select(donor, k = NES),
    
    "concordance" = {m <- d %>%
        select(donor, pathway, NES) %>%
        pivot_wider(names_from = donor, values_from = NES) %>%
        tibble::column_to_rownames("pathway") %>%
        as.matrix()
      
      cons_vec <- apply(m, 1, median, na.rm = TRUE)
      
      data.frame(
        donor = colnames(m),
        k = apply(m, 2, function(v) suppressWarnings(
              cor(v, cons_vec, use = "pairwise.complete.obs", method = "spearman"))))
    },
    
    "ncells" =
      qc %>%
      filter(donor %in% d$donor) %>%
      select(donor, k = n_cells),
    
    "slc7a7" =
      qc %>%
      filter(donor %in% d$donor) %>%
      select(donor, k = mean_g)
  )
  
  col_lvl <- col_key %>%
    arrange(desc(k)) %>%
    pull(donor)
  
  
  ## 셀 값 -------------------------------------------------------------------
  d <- d %>%
    mutate(
      pathway = factor(pathway, levels = row_lvl),
      donor   = factor(donor, levels = col_lvl),
      sig     = padj < SIG_CUT,
      fill_v  = if (GREY_NONSIG)
        ifelse(sig, NES, NA_real_)
      else
        NES
    )
  
  
  ## pathway label -----------------------------------------------------------
  lab_map <- consist %>%
    mutate(
      lab = sprintf(
        "%s  (%d%%)",
        sub("HALLMARK_", "", pathway),
        round(100 * frac_same)
      )
    ) %>%
    select(pathway, lab)
  
  labs_v <- setNames(
    lab_map$lab,
    lab_map$pathway
  )[row_lvl]
  
  
  ## color limit -------------------------------------------------------------
  lim <- max(
    1,
    quantile(abs(d$NES), 0.98, na.rm = TRUE)
  )
  
  
  ## plot --------------------------------------------------------------------
  p <- ggplot(
    d,
    aes(donor, pathway, fill = fill_v)
  ) +
    geom_tile() +
    
    scale_fill_gradient2(
      low = "#2C6FAF",
      mid = "white",
      high = "#C0392B",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      na.value = "grey90",
      name = "NES"
    ) +
    
    scale_x_discrete(
      labels = NULL,
      breaks = NULL
    ) +
    
    scale_y_discrete(
      labels = labs_v
    ) +
    
    labs(
      x = sprintf(
        "donor (n = %d, %s 순 정렬)",
        nlevels(d$donor),
        COL_ORDER
      ),
      y = NULL,
      title = ttl,
      subtitle = sub,
      caption = paste0(
        "각 셀 = donor 내부 세포 수준 GSEA ",
        "(gene ~ SLC7A7 + depth 의 t 통계량 랭킹). ",
        if (GREY_NONSIG)
          "회색 = donor 내 BH FDR ≥ 0.05. "
        else
          "",
        "괄호 = 같은 방향 donor 비율.\n",
        "모든 pathway가 유의하지 않은 donor는 시각화에서 제외함. ",
        "이 그림은 일관성 시각화이며 유의성 검정이 아니다 — ",
        "통계적 주장은 within/between 혼합모형 GSEA(본문)에 있다."
      )
    ) +
    
    theme_paper(base_size = 8) +
    
    theme(
      axis.text.y = element_text(size = 6.4),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.caption = element_text(
        hjust = 0,
        size = 6,
        color = "grey35"
      )
    )
  
  
  if (FACET_BY_DATASET)
    p <- p +
    facet_grid(
      ~ dataset,
      scales = "free_x",
      space = "free_x"
    ) +
    theme(
      strip.text.x = element_text(
        size = 5.5,
        angle = 90
      )
    )
  
  p
}

pA <- make_heatmap(
  "SLC7A7",
  "Donor-wise Hallmark enrichment along within-donor SLC7A7 co-expression",
  sprintf("microglia ≥ %d cells, SLC7A7⁺ ≥ %d cells; depth-adjusted",
          MIN_CELLS_DONOR, MIN_POS_CELLS))
ggsave(file.path(DIR_FIGURES, "Fig6_donor_gsea_heatmap.pdf"),
       pA, width = 8, height = 8.6)
save_emf(pA , "Fig6_donor_gsea_heatmap.emf", legend_p = "right", width = 6, height = 9)

library(devEMF)
emf(file = file.path(DIR_FIGURES, "Fig6_donor_gsea_heatmap.emf"), width = 6, height = 9)

# 그래프 생성 코드
pA

# 파일 저장 완료 (반드시 실행)
dev.off()
pB <- make_heatmap(
  "depth",
  "NEGATIVE CONTROL — same donors, ranked by t(β_depth)",
  "SLC7A7 패널과 유사한 패턴이 보이면 signal 이 depth 교란임을 뜻한다")
ggsave(file.path(DIR_FIGURES, "Fig6_donor_gsea_heatmap_depth.pdf"),
       pB, width = 9.0, height = 8.6)

# (c) 일관성 요약 막대 — 히트맵 옆/아래에 붙일 보조 패널
pC <- consist %>%
  mutate(lab = sub("HALLMARK_", "", pathway),
         signed_frac = ifelse(dir_major == "up", frac_same, -frac_same)) %>%
  ggplot(aes(signed_frac, reorder(lab, signed_frac), fill = dir_major)) +
  geom_col(width = .72) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey60") +
  scale_fill_manual(values = c(up = "#C0392B", down = "#2C6FAF"), guide = "none") +
  scale_x_continuous(labels = function(x) paste0(abs(x) * 100, "%"),
                     limits = c(-1, 1)) +
  labs(x = "같은 방향 donor 비율 (좌 = down, 우 = up)", y = NULL,
       title = "Directional consistency across donors",
       caption = "점선 = 50%(무작위). 이 비율은 기술통계이며 효과크기가 아니다.") +
  theme_paper(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.4),
        plot.caption = element_text(hjust = 0, size = 6, color = "grey35"))
ggsave(file.path(DIR_FIGURES, "Fig6_donor_gsea_consistency.pdf"),
       pC, width = 6.2, height = 8.6)

message("완료: donor별 GSEA 히트맵")

# -----------------------------------------------------------------------------
# 보고 · 해석 가이드
#
# [Methods 에 반드시 들어갈 문장]
#   "Donor-wise GSEA was performed for descriptive visualization only. Donors
#    were included if they contributed ≥100 microglial cells and ≥20 SLC7A7-
#    detected cells (detection rate ≥3%); N donors were excluded (see Table Sx).
#    Within each donor, genes were ranked by the t-statistic of SLC7A7 from a
#    per-gene linear model including sequencing depth as a covariate. Because
#    this analysis is performed at the cell level within donors, it does not
#    constitute an independent statistical test; inferential claims rest on the
#    within/between mixed-effects GSEA."
#
# [그림 캡션에서 피할 표현]
#   × "N/50 pathways were significantly enriched across donors"
#   × "consistent in 92% of donors (p = ...)"   ← sign test 를 주 근거로 승격 금지
#   ○ "the direction of the within-donor effect was shared by 92% of donors"
#
# [읽는 순서]
#   1) Fig6_donor_gsea_qc.csv 에서 탈락 사유 분포 확인
#      — ②(SLC7A7+ 세포 수)가 ①(총 세포 수)보다 많이 탈락시키면,
#        "세포 ≥100" 이라는 기준 자체를 Methods 에서 재서술해야 한다.
#   2) Fig6_donor_gsea_vs_depth.csv 에서 상위 pathway 의 부호 일치 여부 확인
#      — IFN 계열이 depth 패널에서도 양수면 이 그림은 쓸 수 없다.
#   3) FACET_BY_DATASET <- TRUE 로 한 번 더 렌더링
#      — 일관된 donor 들이 특정 dataset 에 몰려 있으면 ICC_dataset = 0.56 이
#        이 그림에도 그대로 들어온 것이므로 캡션에 명시할 것.
#
# [알려진 한계]
#   · donor 간 NES 는 서로 다른 세포 수에서 계산되어 정밀도가 다르다.
#     색의 진하기를 donor 간에 비교하지 말 것(방향의 일치만 읽는다).
#   · 유전자 universe 는 전역 고정이지만 donor 내 0분산 유전자는 제거되므로
#     pathway size 가 donor 마다 미세하게 다르다(long 표의 size 열 참조).
#   · GSEA 를 SLC7A7 로 랭킹한 뒤 SLC7A7 포함 모듈을 보면 순환논리가 되므로
#     universe 에서 self 를 제외했다(2절). Hallmark 에는 SLC7A7 이 없지만
#     GOBP 로 확장할 경우 수송 관련 term 에서 반드시 재확인할 것.
# -----------------------------------------------------------------------------
