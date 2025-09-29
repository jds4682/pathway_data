# 사용자 라이브러리 경로 설정 (쓰기 권한 문제 해결)
user_lib <- Sys.getenv("R_LIBS_USER", unset = "/usr/local/lib/R/site-library")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
.libPaths(user_lib)

# 패키지 설치 함수 정의
install_if_missing <- function(pkg, repo = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    if (repo == "Bioc") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos="https://cloud.r-project.org/")
      }
      BiocManager::install(pkg, update=FALSE, ask=FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
  } else {
    message(paste(pkg, "is already installed."))
  }
}

# 필수 패키지 설치
install_if_missing("clusterProfiler", "Bioc")
install_if_missing("org.Hs.eg.db", "Bioc")
install_if_missing("enrichplot", "Bioc")
install_if_missing("dplyr", "CRAN")
install_if_missing("ggplot2", "CRAN")

print("All required R packages are installed.")
