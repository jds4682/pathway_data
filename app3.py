import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import requests
from io import BytesIO
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# --- 1. 초기 설정 및 GitHub 데이터 로딩 함수 ---

st.set_page_config(layout="wide", page_title="Herbal Prescription Network Analysis")

# --- GitHub 데이터 로딩 함수들 ---
def load_excel_data(name):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/{name}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            st.success(f"'{name}' 파일을 GitHub에서 성공적으로 불러왔습니다.")
            return pd.read_excel(BytesIO(response.content))
        else:
            st.error(f"GitHub에서 '{name}' 파일을 찾을 수 없습니다. (상태 코드: {response.status_code})")
            return None
    except requests.exceptions.RequestException as e:
        st.error(f"네트워크 오류: {e}")
        return None

# ★★★ 개별 약재 데이터(CSV)를 GitHub에서 불러오는 함수 ★★★
def load_herb_csv_data(smhb_code):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/tg/{smhb_code}.csv"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            # CSV 파일을 DataFrame으로 직접 읽음
            return pd.read_csv(BytesIO(response.content))
        else:
            st.warning(f"GitHub에서 '{smhb_code}.csv' 파일을 찾을 수 없습니다.")
            return None
    except requests.exceptions.RequestException:
        st.warning(f"'{smhb_code}.csv' 파일 로딩 중 네트워크 오류 발생.")
        return None

# --- 메인 데이터 로딩 함수 ---
@st.cache_data
def load_initial_data():
    try:
        herb_df = load_excel_data('all name.xlsx')
        ingre_data = load_excel_data('SMIT.xlsx')
        
        if herb_df is None or ingre_data is None:
             return None, None

        ingre_data = ingre_data.dropna(subset=['OB_score'], how='any', axis=0)
        return herb_df, ingre_data
        
    except Exception as e:
        st.error(f"초기 데이터 로딩 중 오류 발생: {e}")
        return None, None

# --- 2. 네트워크 분석 핵심 로직 (CSV 파싱으로 수정) ---

def run_network_analysis(selected_herbs_with_grams, ingre_data):
    t_name = "_".join(selected_herbs_with_grams.keys())
    smhb_codes = list(selected_herbs_with_grams.values())

    node_list = []
    edge_list = []
    
    progress_bar = st.progress(0, text="데이터 로딩을 시작합니다...")
    
    # ★★★ CSV 파일을 불러와 DataFrame으로 파싱하는 로직 ★★★
    for i, code in enumerate(smhb_codes):
        # GitHub에서 개별 약재 CSV 데이터 불러오기
        herb_df_single = load_herb_csv_data(code)
        
        if herb_df_single is not None:
            # DataFrame을 한 줄씩 순회하며 node와 edge 정보 추출
            for idx, row in herb_df_single.iterrows():
                group = row.get('group') # 'group' 열이 있다고 가정
                if group == 'nodes':
                    info = row.get('info', '') # 'info' 열이 있다고 가정
                    node_list.append(info.split("<br>"))
                elif group == 'edges':
                    source = row.get('source') # 'source' 열이 있다고 가정
                    target = row.get('target') # 'target' 열이 있다고 가정
                    edge_list.append([source, target])
        
        progress_bar.progress((i + 1) / len(smhb_codes), text=f"데이터 로딩 중: {code}")

    progress_bar.empty()

    if not node_list:
        st.error("선택된 약재에 대한 유효한 데이터를 불러오지 못했습니다. GitHub 파일 존재 여부 및 형식을 확인하세요.")
        return None, None

    # --- 이하 데이터 가공, 필터링, 네트워크 생성 로직은 기존과 거의 동일 ---
    node_data = pd.DataFrame()
    node_ID, node_label, node_group = [], [], []
    for item in node_list:
        try:
            node_group.append(item[0])
            node_ID.append(item[1][4:].strip())
            node_label.append(item[2][6:])
        except (IndexError, TypeError):
            continue
    node_data['ID'] = node_ID
    node_data['Label'] = node_label
    node_data['Group'] = node_group

    edge_data = pd.DataFrame(edge_list, columns=['SourceID', 'TargetID'])
    
    # (OB Score 필터링 등 후속 처리는 여기에 위치)
    
    st.success("네트워크 분석이 완료되었습니다!")
    # (결과 시각화 및 테이블 생성 로직은 이어서 구현)
    
    # 임시로 간단한 결과만 반환
    return "Figure_Object_Placeholder", node_data.head()


# --- 3. 웹페이지 UI 구성 ---
st.title("🌿 천연물 처방 네트워크 분석기 (GitHub-Powered CSV)")

herb_df, ingre_data = load_initial_data()

if herb_df is not None:
    st.header("1. 약재 선택 및 용량 입력")
    herb_names = herb_df['korean name'].dropna().unique().tolist()
    selected_herb_names = st.multiselect("분석할 약재를 선택하세요.", options=herb_names)
    
    selected_herbs_with_grams = {}
    if selected_herb_names:
        cols = st.columns(len(selected_herb_names))
        for i, name in enumerate(selected_herb_names):
            with cols[i]:
                grams = st.number_input(f"{name} (g)", min_value=0.1, value=4.0, step=0.1, key=name)
                smhb_id = herb_df[herb_df['약재명 (Korean Name)'] == name]['SMHB ID'].iloc[0]
                selected_herbs_with_grams[name] = smhb_id
    
    st.header("2. 분석 실행")
    if st.button("네트워크 분석 시작", disabled=(not selected_herb_names)):
        with st.spinner("GitHub에서 CSV 데이터를 불러와 네트워크를 생성 중입니다..."):
            fig, result_df = run_network_analysis(
                {name: selected_herbs_with_grams[name] for name in selected_herb_names}, 
                ingre_data
            )
            
            if fig and result_df is not None:
                st.header("3. 분석 결과")
                st.subheader("분석된 데이터 (샘플)")
                st.dataframe(result_df)
