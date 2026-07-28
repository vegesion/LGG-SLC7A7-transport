# renv 직접 설정하기 / Set up renv from scratch

이 저장소에는 **임의로 만든 renv.lock 을 넣지 않았습니다.** renv.lock 은 본인 컴퓨터의
실제 설치 상태에서 생성돼야 정확하기 때문입니다. 아래 공식 절차대로 직접 초기화하세요.

## 0. (권장) 먼저 패키지 설치
renv 초기화 시 재다운로드를 줄이려면, 필요한 패키지를 먼저 시스템에 깔아두세요.
```r
source("install_dependencies.R")
```

## 1. renv 초기화 → renv.lock 생성
RStudio에서 `LGG-SLC7A7-microglia.Rproj` 를 연 상태에서:
```r
install.packages("renv")
renv::init()
```
- renv 가 `analysis/`·`R/` 코드의 `library()`/`::` 를 스캔해 필요한 패키지를 파악하고,
  프로젝트 전용 라이브러리(`renv/library/`)로 가져온 뒤 **`renv.lock` 을 생성**합니다.
- 프롬프트가 나오면 기본값(실제 사용하는 패키지만 스냅샷)을 선택하면 됩니다.
- Bioconductor 패키지도 renv 가 자동 인식해 Bioconductor 버전과 함께 기록합니다.

이 시점에 만들어진 `renv.lock` 이 **본인 환경의 정확한 버전 고정본**입니다.

## 2. 평소 사용
```r
renv::status()     # lock 과 현재 라이브러리 차이 확인
renv::snapshot()   # 패키지 추가/업데이트 후 lock 갱신
renv::restore()    # lock 대로 되돌리기 (다른 컴퓨터/협업자와 동기화)
```

## 문제 발생 시
- **`node stack overflow`**: `renv::init()` 의 세션 재시작이 이전에 남은 `.Rprofile`/`renv/`
  와 충돌할 때 납니다. 아래로 완전 초기화 후 다시 시도하세요.
  ```r
  unlink("renv", recursive = TRUE)
  if (file.exists(".Rprofile"))     file.rename(".Rprofile", ".Rprofile.bak")
  if (file.exists(".Rprofile.bak")) file.remove(".Rprofile.bak")
  # → R 재시작(Ctrl+Shift+F10) 후
  renv::init(restart = FALSE)   # restart = FALSE 로 재시작 재귀 자체를 회피
  ```
- 세션이 아예 안 열리면 탐색기에서 `.Rprofile` 이름을 바꾼 뒤 R 을 재시작하세요.
- **다운로드 오류(`readRDS ... unknown input format`)**: 새 세션에서
  `options(download.file.method = "libcurl", timeout = 600)` 후 재시도 (프록시/VPN 확인).

## 참고
renv 는 **선택 사항**입니다. 패키지만 설치돼 있으면 renv 없이도 `analysis/` 파이프라인은
그대로 실행됩니다. renv.lock 은 재현성·협업이 필요할 때 위 절차로 만들면 됩니다.
