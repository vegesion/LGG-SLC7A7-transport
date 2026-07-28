st <- read.csv(file.path(DIR_RESULTS, "Fig3_stratified_cox_idh.csv"))
mut <- subset(st, stratum == "idh=Mutant")
b <- log(mut$HR); se <- (log(mut$upper) - log(mut$lower)) / (2*1.96)
w <- 1/se^2; pooled <- sum(w*b)/sum(w); sep <- sqrt(1/sum(w))
Q <- sum(w*(b-pooled)^2); I2 <- max(0, (Q-(length(b)-1))/Q)*100
c(HR = exp(pooled), lo = exp(pooled-1.96*sep), hi = exp(pooled+1.96*sep),
  p = 2*pnorm(-abs(pooled/sep)), I2 = I2)   # I² 낮으면 "이질성 없이 일관"


"결과 : 메타분석은 훌륭해요. IDH-mutant 통합 HR 1.47 (1.28–1.68), p=3.4e-8, I²=0. I²가 0이라는 건 세 코호트(TCGA_LGG·CGGA_693·CGGA_325)의 효과크기가 통계적으로 구분 안 될 만큼 일치한다는 뜻이에요. 다중 코호트 검증에서 이보다 깔끔하게 나오기 어렵습니다. 이게 Fig 3의 헤드라인이에요."