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
    except requests.exceptions.RequestException as e:
        st.error(f"네트워크 오류: {e}")
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
    Target = [t_name]
    
    progress_bar = st.progress(0, text="약재 데이터를 GitHub에서 로딩 중입니다...")
    smhb_codes = list(selected_herbs_info.values())
    
    for i, code in enumerate(smhb_codes):
        herb_df_single = load_herb_csv_data(code)
        
        # --- ★★★ 오류 수정 부분 ★★★ ---
        # if herb_data: -> if herb_df_single is not None and not herb_df_single.empty:
        if herb_df_single is not None and not herb_df_single.empty:
            # CSV를 JSON 유사 구조(딕셔너리 리스트)로 변환
            herb_json_structure = []
            for _, row in herb_df_single.iterrows():
                group = row.get('group')
                if group == 'nodes':
                    herb_json_structure.append({
                        'group': 'nodes',
                        'data': {'info': row.get('info')}
                    })
                elif group == 'edges':
                    herb_json_structure.append({
                        'group': 'edges',
                        'data': {
                            'source': row.get('source'),
                            'target': row.get('target')
                        }
                    })
            Target.append(herb_json_structure)
        
        progress_bar.progress((i + 1) / len(smhb_codes))
    progress_bar.empty()

    if len(Target) <= 1:
        st.error("선택된 약재에 대한 유효 데이터를 불러오지 못했습니다.")
        return None, None, None, None

    # --- ▼▼▼ 제공해주신 분석 코드 시작 ▼▼▼ ---
    t_name = Target.pop(0)
    node_list = []
    edge_list = []
    for herb in Target:
        for i in range(len(herb)):
            if herb[i]['group'] == 'nodes':
                node_list.append(herb[i]['data']['info'].split("<br>"))
            if herb[i]['group'] == 'edges':
                edge_list.append([herb[i]['data']['source'], herb[i]['data']['target']])

    # data shaping
    node_data = pd.DataFrame()
    node_ID, node_label, node_group = [], [], []
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
    
    edge_data = pd.DataFrame(edge_list, columns=['SourceID', 'TargetID'])

    # OB score 필터링
    st.write("OB Score 기반으로 유효성분 필터링 중...")
    drop_list = []
    valid_ingredients = set(ingre_data['Molecule_name'])
    for i, row in node_data.iterrows():
        if row['Group'] in ['MM symptom', 'TCM symptom']:
            drop_list.append(i)
        elif row['Group'] == 'ingredient' and row['Label'] not in valid_ingredients:
            drop_list.append(i)
    node_data.drop(list(set(drop_list)), axis=0, inplace=True)
    st.write(f"유효성분 필터링 완료. {len(drop_list)}개의 노드 제거.")
    st.warning("주의: 실시간 웹 환경의 제약으로 인해, FDR(q-value) 기반 타겟 필터링은 현재 버전에서 생략되었습니다.")
    
    # 네트워크 생성 및 시각화 (이하 로직은 제공해주신 코드 기반)
    G = nx.Graph()
    for _, row in node_data.iterrows():
        color_map = {'herb': 'orange', 'ingredient': 'green', 'disease': 'yellow', 'target': 'skyblue'}
        G.add_node(row['Label'], ID=row['ID'], Group=row['Group'], color=color_map.get(row['Group'], 'gray'))
    
    disease_label = []
    for _, row in edge_data.iterrows():
        try:
            source_name = node_data.loc[node_data['ID'] == row['SourceID'], 'Label'].iloc[0]
            destination_name = node_data.loc[node_data['ID'] == row['TargetID'], 'Label'].iloc[0]
            G.add_edge(source_name, destination_name)
            if row['TargetID'] and row['TargetID'].startswith('SMDE'):
                disease_label.append(destination_name)
        except IndexError:
            continue

    disease_table = pd.Series(disease_label).value_counts().reset_index()
    disease_table.columns = ['Disease', 'Count']

    node_groups = {g: [n for n, d in G.nodes(data=True) if d.get('Group') == g] for g in ['herb', 'ingredient', 'target', 'disease']}
    
    st.write("네트워크 시각화 생성 중...")
    fig, ax = plt.subplots(figsize=(14, 15))
    shells = [node_groups['herb'], node_groups['ingredient'], node_groups['target'], node_groups['disease']]
    pos = nx.shell_layout(G, shells)
    
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['herb'], node_color='#ff8800', node_size=150, label='Herb', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['ingredient'], node_color='#00d200', node_size=20, label='Ingredient', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['target'], node_color='#ff3367', node_size=20, label='Target', ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_groups['disease'], node_color='#6600ff', node_size=20, label='Disease', ax=ax)
    nx.draw_networkx_edges(G, pos, width=0.2, alpha=0.5, ax=ax)
    
    ax.legend(scatterpoints=1)
    ax.set_title(t_name.replace("_", " + "), fontname='DejaVu Sans', fontsize=16)
    
    return fig, disease_table.head(20), node_data, edge_data

# --- 3. 웹페이지 UI 구성 ---
st.title("🌿 천연물 처방 네트워크 분석기 (GitHub-Powered)")

herb_df, ingre_data = load_initial_data()

if herb_df is not None:
    st.header("1. 약재 선택 및 용량 입력")
    herb_names = herb_df['korean name'].dropna().unique().tolist()
    selected_herb_names = st.multoselect("분석할 약재를 선택하세요.", options=herb_names)
    
    selected_herbs_info = {}
    if selected_herb_names:
        cols = st.columns(len(selected_herb_names))
        for i, name in enumerate(selected_herb_names):
            with cols[i]:
                grams = st.number_input(f"{name} (g)", min_value=0.1, value=4.0, step=0.1, key=name)
                smhb_id = herb_df[herb_df['korean name'] == name]['SMHB_ID'].iloc[0]
                selected_herbs_info[name] = smhb_id
    
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
