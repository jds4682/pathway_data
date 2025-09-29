import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import requests
from io import BytesIO
import pickle

# --- 1. 초기 설정 및 GitHub 데이터 로딩 함수 ---

st.set_page_config(layout="wide", page_title="Herbal Prescription Network Analysis")

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

@st.cache_data
def load_initial_data():
    herb_df = load_excel_data('all name.xlsx')
    ingre_data = load_excel_data('SMIT.xlsx')
    if herb_df is None or ingre_data is None: return None, None
    ingre_data = ingre_data.dropna(subset=['OB_score'], how='any', axis=0)
    return herb_df, ingre_data

# --- 2. 네트워크 분석 핵심 로직 ---
def run_network_analysis(selected_herbs_info, ingre_data):
    
    t_name = "_".join(selected_herbs_info.keys())
    
    # Step 1: 사용자가 선택한 약재의 CSV 파일들을 DataFrame 리스트로 로딩
    Target_DataFrames = []
    progress_bar = st.progress(0, text="약재 데이터를 GitHub에서 로딩 중입니다...")
    smhb_codes = list(selected_herbs_info.values())
    
    for i, code in enumerate(smhb_codes):
        herb_df_single = load_herb_csv_data(code)
        if herb_df_single is not None and not herb_df_single.empty:
            Target_DataFrames.append(herb_df_single)
        progress_bar.progress((i + 1) / len(smhb_codes))
    progress_bar.empty()

    if not Target_DataFrames:
        st.error("선택된 약재에 대한 유효 데이터를 불러오지 못했습니다.")
        return None, None, None, None

    # --- ▼▼▼ 제공해주신 분석 코드의 '의도'를 CSV(DataFrame)에 맞게 구현 ▼▼▼ ---
    
    st.info("로드된 CSV 데이터를 파싱하여 Node와 Edge 목록을 생성합니다...")
    node_list = []
    edge_list = []

    # Step 2: 각 약재의 DataFrame을 순회하며 node_list와 edge_list 생성
    # 이 부분이 CSV 형식에 맞춰 재구성된 핵심 파싱 로직입니다.
    for herb_df in Target_DataFrames:
        for index, row in herb_df.iterrows():
            group = row.get('group')
            if group == 'nodes':
                # info 열의 값이 문자열인지 확인 후 split
                if isinstance(row.get('info'), str):
                    node_list.append(row['info'].split("<br>"))
            elif group == 'edges':
                edge_list.append([row.get('source'), row.get('target')])

    # Step 3: 이하 제공해주신 코드 로직을 기반으로 데이터 가공 및 분석 수행
    # data shaping
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

    # OB score 필터링
    st.info("OB Score 기반으로 유효성분 필터링 중...")
    drop_list = []
    valid_ingredients = set(ingre_data['Molecule_name'])
    for i, row in node_data.iterrows():
        if row['Group'] in ['MM symptom', 'TCM symptom']:
            drop_list.append(i)
        elif row['Group'] == 'ingredient' and row['Label'] not in valid_ingredients:
            drop_list.append(i)
    node_data.drop(list(set(drop_list)), axis=0, inplace=True)
    st.info(f"유효성분 필터링 완료. {len(drop_list)}개의 노드 제거.")
    st.warning("주의: 실시간 웹 환경의 제약으로 인해, FDR(q-value) 기반 타겟 필터링은 현재 버전에서 생략되었습니다.")
    
    # 네트워크 생성 및 시각화
    G = nx.Graph()
    # (이하 네트워크 생성 및 시각화 로직은 이전과 동일)
    for _, row in node_data.iterrows():
        color_map = {'herb': 'orange', 'ingredient': 'green', 'disease': 'yellow', 'target': 'skyblue'}
        G.add_node(row['Label'], ID=row['ID'], Group=row['Group'], color=color_map.get(row['Group'], 'gray'))
    
    disease_label = []
    for _, row in edge_data.iterrows():
        try:
            source_name = node_data.loc[node_data['ID'] == row['SourceID'], 'Label'].iloc[0]
            destination_name = node_data.loc[node_data['ID'] == row['TargetID'], 'Label'].iloc[0]
            G.add_edge(source_name, destination_name)
            if row['TargetID'] and str(row['TargetID']).startswith('SMDE'):
                disease_label.append(destination_name)
        except IndexError:
            continue

    disease_table = pd.Series(disease_label).value_counts().reset_index()
    disease_table.columns = ['Disease', 'Count']

    node_groups = {g: [n for n, d in G.nodes(data=True) if d.get('Group') == g] for g in ['herb', 'ingredient', 'target', 'disease']}
    
    st.info("네트워크 시각화 생성 중...")
    fig, ax = plt.subplots(figsize=(14, 15))
    shells = [node_groups.get('herb', []), node_groups.get('ingredient', []), node_groups.get('target', []), node_groups.get('disease', [])]
    pos = nx.shell_layout(G, shells)
    
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups.get('herb', []), node_color='#ff8800', node_size=150, label='Herb', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups.get('ingredient', []), node_color='#00d200', node_size=20, label='Ingredient', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups.get('target', []), node_color='#ff3367', node_size=20, label='Target', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups.get('disease', []), node_color='#6600ff', node_size=20, label='Disease', ax=ax)
    nx.draw_networkx_edges(G, pos, width=0.2, alpha=0.5, ax=ax)
    
    ax.legend(scatterpoints=1)
    ax.set_title(t_name.replace("_", " + "), fontname='DejaVu Sans', fontsize=16)
    
    return fig, disease_table.head(20), node_data, edge_data


# --- 3. 웹페이지 UI 구성 ---
st.title("🌿 천연물 처방 네트워크 분석기 (GitHub-Powered)")

herb_df, ingre_data = load_initial_data()

if herb_df is not None:
    st.header("1. 약재 선택")
    
    # 'all name.xlsx' 파일의 실제 열 이름을 확인 후 수정 필요
    KOREAN_NAME_COLUMN = 'korean name'
    SMHB_ID_COLUMN = 'SMHB_ID'
    
    try:
        herb_names = herb_df[KOREAN_NAME_COLUMN].dropna().unique().tolist()
        selected_herb_names = st.multiselect("분석할 약재를 선택하세요.", options=herb_names)
        
        selected_herbs_info = {name: herb_df[herb_df[KOREAN_NAME_COLUMN] == name][SMHB_ID_COLUMN].iloc[0] for name in selected_herb_names}
    
        st.header("2. 분석 실행")
        if st.button("네트워크 분석 시작", disabled=(not selected_herb_names)):
            with st.spinner("분석을 실행합니다. 약재 수에 따라 시간이 걸릴 수 있습니다..."):
                fig, disease_df, node_df, edge_df = run_network_analysis(selected_herbs_info, ingre_data)
                
                if fig and disease_df is not None:
                    st.header("3. 분석 결과")
                    st.pyplot(fig)
                    
                    st.subheader("상위 20개 연관 질병")
                    st.dataframe(disease_df)
    
                    st.subheader("결과 데이터 다운로드")
                    col1, col2 = st.columns(2)
                    with col1:
                        st.download_button(
                            label="Node 데이터 다운로드 (CSV)",
                            data=node_df.to_csv(index=False).encode('utf-8-sig'),
                            file_name=f"{'_'.join(selected_herbs_info.keys())}_nodes.csv",
                            mime='text/csv',
                        )
                    with col2:
                        st.download_button(
                            label="Edge 데이터 다운로드 (CSV)",
                            data=edge_df.to_csv(index=False).encode('utf-8-sig'),
                            file_name=f"{'_'.join(selected_herbs_info.keys())}_edges.csv",
                            mime='text/csv',
                        )
    except KeyError as e:
        st.error(f"'{e}' 열을 'all name.xlsx' 파일에서 찾을 수 없습니다. 코드의 열 이름을 확인하세요.")
