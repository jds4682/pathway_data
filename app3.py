import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import json # SMHB 파일을 읽기 위해 추가
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

# --- 1. 초기 설정 및 데이터 로딩 ---

def load_data(name):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/{name}"
    response = requests.get(url)
    print("새로 url 받았음")
    if response.status_code == 200:
        return pd.read_excel(BytesIO(response.content))
        
    else:
        st.error("Could not able to download contents.")
        return None

# 페이지 설정
st.set_page_config(layout="wide", page_title="Herbal Prescription Network Analysis")

# 데이터 파일 경로 (사용자 환경에 맞게 수정)
# 이 파일들은 app.py와 같은 폴더에 있거나, 전체 경로를 지정해야 합니다.
HERB_DB_PATH = r'C:\Users\seoku\Desktop\논문리비전\all name.xlsx' # 약재명 <-> SMHB_ID 매칭 파일
OB_SCORE_PATH = r'C:\Users\seoku\Downloads\SMIT.xlsx' # OB Score 파일 경로
HERB_DATA_DIR = r'C:\tg' # 개별 약재 데이터 폴더 경로

# 데이터 로딩 함수 (캐싱을 사용하여 속도 향상)
@st.cache_data
def load_data():
    try:
        herb_df = load_data('all name.xlsx')
        ingre_data = load_data('SMIT.xlsx')
        ingre_data = ingre_data.dropna(subset=['OB_score'], how='any', axis=0)
        return herb_df, ingre_data
    except FileNotFoundError as e:
        st.error(f"오류: 필수 데이터 파일을 찾을 수 없습니다. 경로를 확인하세요. ({e})")
        return None, None

# --- 2. 네트워크 분석 핵심 로직 (사용자 코드 기반 함수화) ---

def run_network_analysis(selected_herbs_with_grams, ingre_data, herb_data_dir):
    """
    선택된 약재와 용량을 기반으로 네트워크 분석을 수행하는 메인 함수
    """
    t_name = "_".join(selected_herbs_with_grams.keys())
    smhb_codes = list(selected_herbs_with_grams.values())

    # --- 데이터 로딩 및 파싱 ---
    node_list = []
    edge_list = []
    
    # st.write(f"선택된 SMHB 코드: {smhb_codes}")

    for code in smhb_codes:
        file_path = os.path.join(herb_data_dir, f"{code}.json") # 파일 확장자를 .json으로 가정
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                herb_data = json.load(f)
                for item in herb_data:
                    if item.get('group') == 'nodes':
                        node_list.append(item['data']['info'].split("<br>"))
                    elif item.get('group') == 'edges':
                        edge_list.append([item['data']['source'], item['data']['target']])
        except Exception as e:
            st.warning(f"파일을 읽는 중 오류 발생: {file_path} - {e}")
            continue

    # --- 데이터 가공 (Shaping) ---
    node_data = pd.DataFrame()
    edge_data = pd.DataFrame()
    node_ID = []
    node_label = []
    node_group = []
    for i in range(len(node_list)):
        try:
            node_ID.append(node_list[i][1][4:].strip())
            node_label.append(node_list[i][2][6:])
            node_group.append(node_list[i][0])
        except IndexError:
            continue
    node_data['ID'] = node_ID
    node_data['Label'] = node_label
    node_data['Group'] = node_group

    edge_data['SourceID'] = [edge[0] for edge in edge_list]
    edge_data['TargetID'] = [edge[1] for edge in edge_list]

    # --- 필터링 (OB Score, FDR) ---
    # (사용자 코드의 필터링 로직을 여기에 통합합니다)
    # OB Score 필터링
    valid_ingredients = set(ingre_data['Molecule_name'])
    drop_indices = []
    for i, row in node_data.iterrows():
        if row['Group'] == 'ingredient' and row['Label'] not in valid_ingredients:
            drop_indices.append(i)
    node_data.drop(drop_indices, axis=0, inplace=True)

    # (FDR 필터링 로직은 매우 복잡하고 모든 파일에 대한 사전 계산이 필요하므로,
    # 여기서는 단순화하거나, 미리 계산된 결과를 사용하는 방식으로 변경하는 것을 권장합니다.)
    # 지금은 이 부분을 생략하고 진행합니다.

    # --- 네트워크 생성 및 시각화 ---
    G = nx.Graph()
    for i, row in node_data.iterrows():
        color_map = {'herb': 'orange', 'ingredient': 'green', 'disease': 'yellow', 'target': 'skyblue'}
        G.add_node(row['Label'], ID=row['ID'], Group=row['Group'], color=color_map.get(row['Group'], 'gray'))

    for i, row in edge_data.iterrows():
        try:
            source_name = node_data[node_data['ID'] == row['SourceID']]['Label'].iloc[0]
            destination_name = node_data[node_data['ID'] == row['TargetID']]['Label'].iloc[0]
            G.add_edge(source_name, destination_name)
        except IndexError:
            continue

    # --- 질병 카운트 ---
    disease_label = edge_data[edge_data['TargetID'].str.startswith('SMDE', na=False)]['TargetID'].map(
        node_data.set_index('ID')['Label']
    ).dropna()
    disease_table = disease_label.value_counts().reset_index()
    disease_table.columns = ['Disease', 'Count']

    # --- 시각화 ---
    fig, ax = plt.subplots(figsize=(14, 15))
    node_groups = {g: [n for n, d in G.nodes(data=True) if d.get('Group') == g] for g in ['herb', 'ingredient', 'target', 'disease']}
    shells = [node_groups['herb'], node_groups['ingredient'], node_groups['target'], node_groups['disease']]
    pos = nx.shell_layout(G, shells)
    
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['herb'], node_color='#ff8800', node_size=100, label='Herb', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['ingredient'], node_color='#00d200', node_size=10, label='Ingredient', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['target'], node_color='#ff3367', node_size=10, label='Target', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['disease'], node_color='#6600ff', node_size=10, label='Disease', ax=ax)
    nx.draw_networkx_edges(G, pos, width=0.1, ax=ax)
    
    ax.legend(scatterpoints=1)
    ax.set_title(t_name, fontname='DejaVu Sans', fontsize=16)
    
    return fig, disease_table.head(20) # 상위 20개 질병만 반환

# --- 3. 웹페이지 UI 구성 ---

st.title("🌿 천연물 처방 네트워크 분석기")

# 데이터 로드
herb_df, ingre_data = load_data()

if herb_df is not None:
    # 약재 선택
    st.header("1. 약재 선택")
    herb_names = herb_df['약재명 (Korean Name)'].dropna().unique().tolist()
    selected_herb_names = st.multiselect("분석할 약재를 선택하세요.", options=herb_names)
    
    # 용량 입력
    st.header("2. 용량 입력 (그램)")
    selected_herbs_with_grams = {}
    if selected_herb_names:
        cols = st.columns(len(selected_herb_names))
        for i, name in enumerate(selected_herb_names):
            with cols[i]:
                grams = st.number_input(f"{name} (g)", min_value=0.1, value=4.0, step=0.1)
                smhb_id = herb_df[herb_df['약재명 (Korean Name)'] == name]['SMHB ID'].iloc[0]
                selected_herbs_with_grams[name] = smhb_id
    else:
        st.info("먼저 분석할 약재를 선택해주세요.")

    # 분석 실행
    st.header("3. 분석 실행")
    if st.button("네트워크 분석 시작", disabled=(not selected_herb_names)):
        with st.spinner("네트워크를 생성하고 분석 중입니다. 잠시만 기다려주세요..."):
            # 분석 함수 실행
            fig, disease_df = run_network_analysis(
                {k: selected_herbs_with_grams[k] for k in selected_herb_names}, 
                ingre_data, 
                HERB_DATA_DIR
            )
            
            # 결과 출력
            st.header("4. 분석 결과")
            st.subheader("네트워크 시각화")
            st.pyplot(fig)
            
            st.subheader("상위 20개 연관 질병")
            st.dataframe(disease_df)
