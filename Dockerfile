# 1. 베이스 이미지 설정 (Python 3.9 환경에서 시작)
FROM python:3.9-slim

# 2. 작업 디렉토리 설정
WORKDIR /app

# 3. 시스템 패키지 설치 (R, libcurl 등 packages.txt의 역할)
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# --- ★★★ 수정된 부분 ★★★ ---
# 4. R 라이브러리 경로를 영구적인 환경 변수로 설정
# 이 환경 변수는 install_packages.R과 app.py의 R 세션 모두에서 사용됩니다.
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# 5. Python 라이브러리 설치
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 6. R 라이브러리 설치용 스크립트 복사 및 실행
COPY install_packages.R ./
RUN Rscript install_packages.R

# 7. 나머지 앱 소스 코드 전체 복사
COPY . .

# 8. 앱 실행 명령어
EXPOSE 8501
# CMD 명령어의 파일 이름을 'app.py'로 수정
CMD ["streamlit", "run", "app3.py"]
