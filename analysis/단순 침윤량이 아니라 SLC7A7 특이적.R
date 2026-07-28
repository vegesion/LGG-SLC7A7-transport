"SLC7A7이 일반적인 myeloid 마커들보다 더 나은가?"
"단순 침윤량이 아니라 SLC7A7 특이적"



mk <- intersect(c("CD68","AIF1","CSF1R","C1QA","C1QB","TYROBP","ITGAM","PTPRC"), rownames(expr))
mye_score <- colMeans(expr[mk, , drop = FALSE])

# (a) SLC7A7 이 generic myeloid 함량을 넘어 추가 정보를 주는가
d <- data.frame(time = clin$time, event = clin$event,
                gene = scale(expr[GENE_OF_INTEREST, ])[,1], mye = scale(mye_score)[,1],
                age = clin$age, grade = clin$grade, idh = clin$idh)
summary(coxph(Surv(time, event) ~ gene + mye + age + grade + idh, data = d))


"Call:
coxph(formula = Surv(time, event) ~ gene + mye + age + grade + 
    idh, data = d)

  n= 377, number of events= 52 
   (결측으로 인하여 139개의 관측치가 삭제되었습니다.)

            coef exp(coef) se(coef)      z Pr(>|z|)    
gene     0.86317   2.37066  0.37613  2.295  0.02174 *  
mye     -0.28637   0.75099  0.36309 -0.789  0.43030    
age      0.06645   1.06871  0.01337  4.971 6.66e-07 ***
gradeG3  0.49248   1.63637  0.32591  1.511  0.13076    
idhWT    1.16171   3.19539  0.38375  3.027  0.00247 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

        exp(coef) exp(-coef) lower .95 upper .95
gene        2.371     0.4218    1.1342     4.955
mye         0.751     1.3316    0.3686     1.530
age         1.069     0.9357    1.0411     1.097
gradeG3     1.636     0.6111    0.8639     3.100
idhWT       3.195     0.3130    1.5062     6.779

Concordance= 0.863  (se = 0.027 )
Likelihood ratio test= 84.79  on 5 df,   p=<2e-16
Wald test            = 74.71  on 5 df,   p=1e-14
Score (logrank) test = 112.9  on 5 df,   p=<2e-16"


# (b) 각 myeloid 마커와 head-to-head (HR 비교)
sapply(c(GENE_OF_INTEREST, mk), function(g) {
  x <- scale(expr[g, ])[,1]
  s <- summary(coxph(Surv(clin$time, clin$event) ~ x + clin$age + clin$idh))
  c(HR = s$conf.int[1,1], p = s$coefficients[1,5])
})

"         SLC7A7      CD68       AIF1      CSF1R         C1QA        C1QB      TYROBP       ITGAM        PTPRC
HR 1.731035e+00 1.1559129 1.33425416 1.32889526 1.5499057288 1.608105830 1.532701576 1.473010471 1.5309255657
p  3.139045e-05 0.2488077 0.02194851 0.01624726 0.0003784196 0.000142667 0.001082123 0.002326629 0.0003148134"



# CGGA에서도 맞나 검증 -----------------------------------------------------------


for (nm in c("CGGA_693", "CGGA_325")) {
  co <- load_cohort(nm); cl <- co$clin
  mk <- intersect(c("SLC7A7","CD68","AIF1","CSF1R","C1QA","C1QB","TYROBP","ITGAM","PTPRC"),
                  rownames(co$expr))
  r <- sapply(mk, function(g) {
    s <- summary(coxph(Surv(cl$time, cl$event) ~ scale(co$expr[g,])[,1] + cl$age + cl$idh))
    c(HR = s$conf.int[1,1], p = s$coefficients[1,5]) })
  cat("\n==", nm, "==\n"); print(round(t(r), 4))
}



"[cohort] CGGA_693
  ID 매칭: CGGA_ID (exact) → 693 샘플
  샘플 693 | 생존정보 유효 657 | IDH 642

== CGGA_693 ==
           HR      p
SLC7A7 1.2251 0.0001
CD68   1.0715 0.1878
AIF1   1.0539 0.2945
CSF1R  1.0198 0.7028
C1QA   1.1431 0.0067
C1QB   1.1546 0.0040
TYROBP 1.1165 0.0282
ITGAM  1.0727 0.2122
PTPRC  1.2176 0.0004
[cohort] CGGA_325
  ID 매칭: CGGA_ID (exact) → 325 샘플
  샘플 325 | 생존정보 유효 313 | IDH 324

== CGGA_325 ==
           HR      p
SLC7A7 1.6025 0.0000
CD68   1.1511 0.0431
AIF1   1.3831 0.0000
CSF1R  1.1804 0.0186
C1QA   1.5549 0.0000
C1QB   1.5379 0.0000
TYROBP 1.5425 0.0000
ITGAM  1.3260 0.0005
PTPRC  1.3692 0.0000"