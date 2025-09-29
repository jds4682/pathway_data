import streamlit as st
import pandas as pd
import requests
from io import BytesIO

# --- 1. 초기 설정 및 GitHub 데이터 로딩 함수 ---

st.set_page_config(layout="wide", page_title="GSEA Pre-processing Tool")

# GitHub에서 데이터(엑셀, CSV)를 불러오는 함수들
@st.cache_data
def load_excel_data(name):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/{name}"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            return pd.read_excel(BytesIO(response.content))
        else:
            st.error(f"GitHub에서 '{name}' 파일을 찾을 수 없습니다.")
            return None
    except Exception as e:
        st.error(f"'{name}' 파일 로딩 중 오류 발생: {e}")
        return None

@st.cache_data
def load_herb_csv_data(smhb_code):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/{smhb_code}.csv"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200:
            return pd.read_csv(BytesIO(response.content))
        else:
            st.warning(f"GitHub에서 '{smhb_code}.csv' 파일을 찾을 수 없습니다.")
            return None
    except requests.exceptions.RequestException:
        st.warning(f"'{smhb_code}.csv' 파일 로딩 중 네트워크 오류 발생.")
        return None

# 앱 시작 시 필요한 약재 목록 데이터를 로드하는 함수
@st.cache_data
def load_initial_data():
    herb_df = load_excel_data('all name.xlsx')
    return herb_df

# --- 2. GSEA 전처리 핵심 로직 (사용자 코드 통합) ---

def process_for_gsea(prescription_name, selected_herbs_info, herb_weights):
    
    data_list = []
    
    progress_bar = st.progress(0, text="약재 데이터를 GitHub에서 로딩 및 처리 중입니다...")
    
    # selected_herbs_info는 {'한글약재명': 'SMHB코드'} 형태의 딕셔너리
    # herb_weights는 {'한글약재명': 그람수} 형태의 딕셔너리
    smhb_codes = list(selected_herbs_info.values())
    
    # --- ▼▼▼ 제공해주신 분석 코드의 의도를 반영한 로직 ▼▼▼ ---
    for i, (herb_name, herb_code) in enumerate(selected_herbs_info.items()):
        
        # GitHub에서 개별 약재 CSV 데이터 로드
        df = load_herb_csv_data(herb_code)
        
        if df is None or df.empty:
            st.warning(f"'{herb_name}'({herb_code}) 데이터를 찾을 수 없거나 비어있어 건너뜁니다.")
            continue

        # 숫자가 아닌 값이 있는지 확인하고 제외
        df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
        df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]

        # 가중치 가져오기 (사용자가 입력한 그람 수)
        weight = herb_weights.get(herb_name, 1.0)

        # 새로운 데이터프레임에 저장
        for _, row in df.iterrows():
            gene = row['Gene symbol']
            value = float(row['Value'])
            score = value * weight  # 가중치(그람 수) 적용
            data_list.append([herb_name, gene, score]) # herb_name을 사용해 가독성 높임
            
        progress_bar.progress((i + 1) / len(smhb_codes))
    
    progress_bar.empty()
    
    if not data_list:
        st.error("처리할 데이터가 없습니다.")
        return None

    # 최종 데이터프레임 생성
    output_df = pd.DataFrame(data_list, columns=['Herb', 'GeneSymbol', 'Score'])

    return output_df

# --- 3. 웹페이지 UI 구성 ---
st.title("🌿 GSEA 분석용 전처리 파일 생성기")
st.info("이 앱은 R에서 GSEA 분석을 수행하기 전에 필요한 `_processed.csv` 파일을 생성합니다.")

herb_df = load_initial_data()

if herb_df is not None:
    st.header("1. 처방 구성하기")
    
    # all name.xlsx 파일의 실제 열 이름을 확인 후 수정 필요
    KOREAN_NAME_COLUMN = 'korean name'
    SMHB_ID_COLUMN = 'SMHB_ID'
    
    try:
        herb_names = herb_df[KOREAN_NAME_COLUMN].dropna().unique().tolist()
        
        # 약재 선택
        selected_herb_names = st.multiselect("분석에 포함할 약재를 선택하세요.", options=herb_names)
        
        selected_herbs_info = {}
        herb_weights = {}

        if selected_herb_names:
            st.subheader("용량 입력 (단위: g)")
            # 3열 레이아웃으로 용량 입력 UI 구성
            num_columns = 3
            cols = st.columns(num_columns)
            
            for i, name in enumerate(selected_herb_names):
                with cols[i % num_columns]:
                    grams = st.number_input(f"{name}", min_value=0.1, value=4.0, step=0.1, key=name)
                    
                    smhb_id = herb_df[herb_df[KOREAN_NAME_COLUMN] == name][SMHB_ID_COLUMN].iloc[0]
                    selected_herbs_info[name] = smhb_id
                    herb_weights[name] = grams
            
            st.divider()
            st.header("2. 분석 및 파일 다운로드")
            
            # 처방 이름 입력
            prescription_name_input = st.text_input("저장할 처방의 영문 이름을 입력하세요 (예: My_Prescription):")
            
            if st.button("전처리 파일 생성 시작", disabled=(not prescription_name_input)):
                with st.spinner("파일을 생성 중입니다..."):
                    
                    # GSEA 전처리 함수 실행
                    result_df = process_for_gsea(selected_herbs_info, herb_weights)
                    
                    if result_df is not None:
                        st.success("파일 생성이 완료되었습니다! 아래에서 결과를 확인하고 다운로드하세요.")
                        
                        st.subheader("생성된 데이터 미리보기 (상위 10개)")
                        st.dataframe(result_df.head(10))
                        
                        # 다운로드 버튼
                        st.download_button(
                            label=f"📥 {prescription_name_input}_processed.csv 다운로드",
                            data=result_df.to_csv(index=False).encode('utf-8-sig'),
                            file_name=f"{prescription_name_input}_processed.csv",
                            mime='text/csv',
                        )

    except KeyError as e:
        st.error(f"'{e}' 열을 'all name.xlsx' 파일에서 찾을 수 없습니다. 코드의 열 이름을 확인하세요.")
