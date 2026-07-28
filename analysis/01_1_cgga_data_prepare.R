


for (proj in c("CGGA_693", "CGGA_325")) {

  out <- read.table(file.path(DIR_DATA_PROC, paste0(proj, "_rawcounts.txt")))   # 이미 있으면 재다운로드 안 함
  
  colnames(out)<- out[1,]
  out <- out[-1,]
  rownames(out) <- out[,1]
  out <- out[,-1]
  out[] <- lapply(out, as.numeric)
  
  norm_csv <- file.path(DIR_DATA_PROC, paste0(proj, "_norm.txt"))
  if (!file.exists(norm_csv)) {
    norm <- normalize_logcpm(out, prior_count = 1)
    write.csv(norm, norm_csv, row.names = FALSE)
    message("[saved] ", norm_csv)
  } else {
    message("[cache] norm csv 이미 있음 → ", norm_csv)
  }
}

message("완료: CGGA 준비. 재실행 시 로컬 저장본을 사용하므로 재다운로드하지 않습니다.")

