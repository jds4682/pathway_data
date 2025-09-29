# 1. 일부러 에러를 발생시켜 Dockerfile이 실행되는지 확인하는 테스트 명령어
RUN echo "This is a test to see if the Dockerfile is being executed." && exit 1

# --- 이하 내용은 테스트 중에는 실행되지 않습니다 ---
FROM python:3.9-slim
WORKDIR /app
# ... (기존의 나머지 Dockerfile 내용) ...
