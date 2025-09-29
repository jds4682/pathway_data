# 1. 베이스 이미지 설정 (Python 3.9 환경에서 시작)
FROM python:3.9-slim

# 2. 작업 디렉토리 설정
WORKDIR /app

# 3. 시스템 패키지 설치 (R, libcurl 등 packages.txt의 역할)
# noninteractive 설정을 통해 설치 중 묻는 질문에 자동으로 'yes'로 답함
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# 4. Python 라이브러리 설치
# requirements.txt 파일을 먼저 복사하여, 이 파일이 변경될 때만 라이브러리를 다시 설치하도록 함 (캐싱 효율)
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 5. R 라이브러리 설치용 스크립트 복사 및 실행
# 이 단계에서 시간이 가장 오래 걸리지만, 배포 시 단 한 번만 실행됨
COPY install_packages.R ./
RUN Rscript install_packages.R

# 6. 나머지 앱 소스 코드 전체 복사
COPY . .

# 7. 앱 실행 명령어
# EXPOSE: 8501 포트를 외부에 노출
# CMD: 컨테이너가 시작될 때 실행할 명령어
EXPOSE 8501
CMD ["streamlit", "run", "app3.py"]
