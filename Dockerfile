# 1. 베이스 이미지 설정 (Python 3.9 환경에서 시작)
FROM python:3.9-slim

# 2. 작업 디렉토리 설정
WORKDIR /app

# 3. 시스템 패키지 설치 (R 및 데이터 과학용 필수 라이브러리)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    g++ \
    make \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libfontconfig1-dev \
    libsqlite3-dev \
    unixodbc-dev \
    && rm -rf /var/lib/apt/lists/*

# --- ★★★ 수정된 부분 ★★★ ---
# 4. R 설치 경로를 영구적인 환경 변수로 설정
#    이것이 pip가 rpy2를 빌드할 때 R을 찾을 수 있게 해주는 핵심입니다.
ENV R_HOME=/usr/lib/R
# R 라이브러리 경로도 함께 설정
ENV R_LIBS_USER=/app/r_libs

# 5. R 패키지를 설치할 폴더를 미리 생성하고 권한 부여
RUN mkdir -p ${R_LIBS_USER} && chmod -R 777 ${R_LIBS_USER}

# 6. Python 라이브러리 설치
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 7. R 라이브러리 설치
COPY install_packages.R ./
RUN Rscript install_packages.R

# 8. 나머지 앱 소스 코드 전체 복사
COPY . .

# 9. 앱 실행 명령어
EXPOSE 8501
CMD ["streamlit", "run", "app.py"]

