# TCGA 다운로드 안정화 / Stabilizing GDC downloads

## 증상별 정리
- **`chunks download was not correct` 무한 반복 + "Downloading: 11 MB → completed"**
  → 다운로드가 중간에 **잘림(truncation)**. 청크 실제 크기보다 훨씬 작은 지점에서 끊김.
  원인은 대개 **백신(안티바이러스)·방화벽·프록시**가 tar.gz 다운로드를 검사하며 연결을 끊음.
- **오래 걸리다 실패** → 타임아웃. `options(timeout = 10000)` 로 해결(코드에 반영됨).

## 방법 1: api 방식 + 아주 작은 청크 (빠른 우회)
파일 1개가 ~4MB 이므로 청크를 1~2개로 줄이면 각 청크가 컷오프 아래라 통과할 수 있습니다.
`config/config.R`:
```r
TCGA_DOWNLOAD_METHOD <- "api"
TCGA_FILES_PER_CHUNK <- 1     # 1개씩 (가장 안전)
```
그리고:
```r
options(timeout = 10000, download.file.method = "libcurl")
unlink(DIR_GDC, recursive = TRUE)     # 깨진 부분 삭제 후
source(here::here("analysis", "01_tcga_download_prepare.R"))
```
+ 가능하면 **백신 실시간검사 잠시 끄기**, **다른 네트워크(휴대폰 핫스팟)** 시도.

## 방법 2: gdc-client (근본 해결, 권장)
GDC 공식 다운로드 도구. 재개 가능하고 대용량에 안정적이라 트렁케이션 문제에 강합니다.

1. 내려받기: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
   (Windows 용 zip 다운로드 → 압축 해제)
2. `gdc-client.exe` 를 PATH 에 추가하거나, R 에서 해당 폴더를 PATH 에 등록:
   ```r
   old <- Sys.getenv("PATH")
   Sys.setenv(PATH = paste(old, "C:/tools/gdc-client", sep = ";"))  # exe 있는 폴더
   system("gdc-client --version")   # 인식되는지 확인
   ```
3. `config/config.R` 에서 방식 변경:
   ```r
   TCGA_DOWNLOAD_METHOD <- "client"
   ```
4. 실행:
   ```r
   unlink(DIR_GDC, recursive = TRUE)
   source(here::here("analysis", "01_tcga_download_prepare.R"))
   ```

## 참고
- 어떤 방식이든 **한 번 성공하면** 결과가 `data/processed/*_counts_clinical.rds` 로 저장되고,
  이후에는 재다운로드하지 않습니다(캐시). 다운로드 일시는 `*_download_info.txt` 에 기록됩니다.
- 방법 1로 받다가 자꾸 특정 청크에서만 막히면, 그 세션을 끊지 말고 그대로 두면
  TCGAbiolinks 가 같은 청크를 재시도합니다. 그래도 안 되면 방법 2로 전환하세요.
