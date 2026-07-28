# =============================================================================
# R/cohorts.R  ---  다중 코호트 로더 / multi-cohort loading & harmonisation
# -----------------------------------------------------------------------------
# TCGA(LGG/GBM), CGGA(693/325), GEO 를 동일 구조로 표준화:
#   list(name, expr = genes×samples(log), clin = data.frame(sample,time,event,age,
#        grade,idh,codel,mgmt,cohort))
# ★ 샘플 ID 는 코호트마다 후보가 여럿(barcode/sample/rownames…)이라, 발현행렬
#   컬럼명과 실제로 겹치는 것을 align_expr_clin() 이 자동 선택한다.
# =============================================================================

# ── 발현행렬 로드 (csv/txt 자동, 필요시 log2 변환) ──────────────────────────
load_cohort_expr <- function(file, log_transform = NA) {
  path <- require_data(file)
  ext  <- tolower(tools::file_ext(path))
  m <- if (ext %in% c("txt", "tsv")) {
    utils::read.csv(path, check.names = FALSE, row.names = 1)
  } else {
    d <- utils::read.csv(path, check.names = FALSE); rownames(d) <- d[[1]]; d[, -1, drop = FALSE]
  }
  m <- as.matrix(m); mode(m) <- "numeric"
  if (is.na(log_transform)) log_transform <- max(m, na.rm = TRUE) > 50
  if (log_transform) { message("  log2(x+1) 변환 적용"); m <- log2(m + 1) }
  m[!duplicated(rownames(m)), , drop = FALSE]
}

# ── 임상표 표준화 ───────────────────────────────────────────────────────────
# 여기서는 ID 를 확정하지 않고 후보(.id_*)를 모두 보존한다.
harmonise_clinical <- function(clin, cohort) {
  n <- nrow(clin)
  g <- function(...) { for (nm in c(...)) if (nm %in% names(clin)) return(clin[[nm]]); rep(NA, n) }

  time  <- suppressWarnings(as.numeric(g("paper_Survival..months.", "OS", "os_months",
                                         "Survival_months", "overall_survival_months")))
  ev_raw <- g("paper_Vital.status..1.dead.", "Censor (alive=0; dead=1)", "status",
              "vital_status_num", "Censor", "event")
  event <- suppressWarnings(as.numeric(ev_raw))
  if (all(is.na(event)) && !all(is.na(ev_raw)))
    event <- as.integer(grepl("dead", as.character(ev_raw), ignore.case = TRUE))

  # ★ fallback: paper_* 가 없는 TCGA → vital_status + days_to_death/follow_up
  if (all(is.na(time)) || all(is.na(event))) {
    vs  <- g("vital_status")
    d2d <- suppressWarnings(as.numeric(g("days_to_death")))
    d2f <- suppressWarnings(as.numeric(g("days_to_last_follow_up", "days_to_last_followup")))
    if (!all(is.na(vs))) {
      dead <- grepl("dead", as.character(vs), ignore.case = TRUE)
      time  <- ifelse(dead, d2d, d2f) / 30.44
      event <- as.integer(dead)
      message("  생존정보: vital_status/days_to_* 로부터 구성")
    }
  }
  if (grepl("^CGGA", cohort) && !all(is.na(time)) && stats::median(time, na.rm = TRUE) > 200)
    time <- time / 30.44

  out <- data.frame(
    time = time, event = event,
    age   = suppressWarnings(as.numeric(g("age_at_diagnosis", "Age", "age", "age_at_index"))),
    grade = as.character(g("tumor_grade", "Grade", "paper_Grade", "grade", "who_grade")),
    idh   = as.character(g("paper_IDH.status", "IDH_mutation_status", "IDH.status", "idh")),
    codel = as.character(g("paper_X1p.19q.codeletion", "1p19q_codeletion_status", "codel")),
    mgmt  = as.character(g("paper_MGMT.promoter.status", "MGMTp_methylation_status", "mgmt")),
    cohort = cohort, stringsAsFactors = FALSE)

  out$idh <- dplyr::case_when(grepl("mut", out$idh, ignore.case = TRUE) ~ "Mutant",
                              grepl("wt|wild", out$idh, ignore.case = TRUE) ~ "WT", TRUE ~ NA_character_)
  out$codel <- dplyr::case_when(grepl("non|no", out$codel, ignore.case = TRUE) ~ "non-codel",
                                grepl("codel|yes", out$codel, ignore.case = TRUE) ~ "codel", TRUE ~ NA_character_)
  out$mgmt <- dplyr::case_when(grepl("un", out$mgmt, ignore.case = TRUE) ~ "unmethylated",
                               grepl("meth", out$mgmt, ignore.case = TRUE) ~ "methylated", TRUE ~ NA_character_)
  if (!all(is.na(out$age)) && stats::median(out$age, na.rm = TRUE) > 200) out$age <- out$age / 365.25
  out$grade <- gsub("^(WHO ?)?", "", toupper(out$grade))

  # ID 후보 보존
  for (nm in intersect(c("barcode","sample","sample_id","submitter_id","Sample",
                         "CGGA_ID","sampleID","patient","PatientID"), names(clin)))
    out[[paste0(".id_", nm)]] <- as.character(clin[[nm]])
  out$.id_rownames <- rownames(clin)
  out
}

# ── 발현행렬 ↔ 임상표 정렬 (ID 자동 선택) ───────────────────────────────────
align_expr_clin <- function(expr, clin, cohort = "") {
  target  <- colnames(expr)
  id_cols <- grep("^\\.id_|CGGA_ID", names(clin), value = TRUE)
  if (!length(id_cols)) stop("[", cohort, "] 임상표에 샘플 ID 후보가 없습니다.", call. = FALSE)

  best <- list(n = -1)
  for (ic in id_cols) {
    ids <- clin[[ic]]
    trials <- list(list(t = target, note = "exact"))
    for (k in c(12, 15, 16))                                  # TCGA aliquot → sample/patient
      trials[[length(trials) + 1]] <- list(t = substr(target, 1, k), note = paste0("trunc", k))
    for (tr in trials) {
      n_ov <- length(intersect(unique(ids), unique(tr$t)))
      if (n_ov > best$n) best <- list(n = n_ov, col = ic, target = tr$t, ids = ids, note = tr$note)
    }
  }
  if (best$n < 1) stop(sprintf(
    "[%s] 임상표와 발현행렬의 샘플 ID 가 일치하지 않습니다.\n  발현 예: %s\n  임상 예: %s",
    cohort, paste(utils::head(target, 2), collapse = ", "),
    paste(utils::head(clin[[id_cols[1]]], 2), collapse = ", ")), call. = FALSE)
  message(sprintf("  ID 매칭: %s (%s) → %d 샘플",
                  sub("^\\.id_", "", best$col), best$note, best$n))

  keep <- best$target %in% best$ids & !duplicated(best$target)
  expr2 <- expr[, keep, drop = FALSE]
  clin2 <- clin[match(best$target[keep], best$ids), , drop = FALSE]
  clin2$sample <- colnames(expr2); rownames(clin2) <- NULL
  clin2 <- clin2[, c("sample", setdiff(names(clin2), c("sample", id_cols))), drop = FALSE]
  list(expr = expr2, clin = clin2)
}

# ── 코호트 1개 로드 ─────────────────────────────────────────────────────────
load_cohort <- function(name, cohorts = COHORTS) {
  spec <- cohorts[[name]]
  if (is.null(spec)) stop("정의되지 않은 코호트: ", name, call. = FALSE)
  message("[cohort] ", name)
  if (grepl("^TCGA", name)) {
    tag <- sub("TCGA_", "", name)
    clin_raw <- readRDS(file.path(DIR_DATA_PROC, paste0(tag, "_counts_clinical.rds")))$clinical
    expr <- load_cohort_expr(spec$expr, log_transform = FALSE)
  } else {
    expr <- load_cohort_expr(spec$expr)
    clin_raw <- utils::read.delim(require_data(spec$clin), check.names = FALSE)
  }
  clin <- harmonise_clinical(clin_raw, name)
  al   <- align_expr_clin(expr, clin, name)

  n_ok <- sum(stats::complete.cases(al$clin[, c("time","event")]) & al$clin$time > 0, na.rm = TRUE)
  message(sprintf("  샘플 %d | 생존정보 유효 %d | IDH %d",
                  ncol(al$expr), n_ok, sum(!is.na(al$clin$idh))))
  if (n_ok < 10) warning("[", name, "] 생존정보가 거의 없습니다. 임상 컬럼명을 확인하세요.", call. = FALSE)
  c(list(name = name), al)
}

load_all_cohorts <- function(cohort_names = names(COHORTS)) {
  out <- lapply(cohort_names, function(n)
    tryCatch(
      load_cohort(n),
      error = function(e) {
        message("  건너뜀(", n, "): ", conditionMessage(e))
        NULL
      }
    )
  )
  
  names(out) <- cohort_names
  out[!sapply(out, is.null)]
}

cohort_gene <- function(co, gene = GENE_OF_INTEREST) {
  if (!gene %in% rownames(co$expr)) return(rep(NA_real_, ncol(co$expr)))
  as.numeric(co$expr[gene, ])
}
