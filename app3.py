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
    
    # --- ★★★ 디버깅 코드 추가 시작 ★★★ ---
    st.subheader("🛠️ 데이터 로딩 상태 (디버깅)")
    # 로그 메시지를 담을 컨테이너 생성
    log_container = st.container(border=True)
    # --- ★★★ 디버깅 코드 추가 끝 ★★★ ---

    # Step 1: 사용자가 선택한 약재의 CSV 파일들을 DataFrame 리스트로 로딩
    Target_DataFrames = []
    progress_bar = st.progress(0, text="약재 데이터를 GitHub에서 로딩 중입니다...")
    smhb_codes = list(selected_herbs_info.values())
    
    for i, code in enumerate(smhb_codes):
        # --- ★★★ 디버깅 코드 추가 시작 ★★★ ---
        log_container.write(f"🔄 `{code}` 파일 로딩 시도 중...")
        # --- ★★★ 디버깅 코드 추가 끝 ★★★ ---

        herb_df_single = load_herb_csv_data(code)
        
        if herb_df_single is not None and not herb_df_single.empty:
            # --- ★★★ 디버깅 코드 추가 시작 ★★★ ---
            log_container.success(f"✅ 성공: `{code}` 파일 로드 완료. (총 {len(herb_df_single)} 행)")
            # --- ★★★ 디버깅 코드 추가 끝 ★★★ ---
            Target_DataFrames.append(herb_df_single)
        else:
            # --- ★★★ 디버깅 코드 추가 시작 ★★★ ---
            log_container.error(f"❌ 실패: `{code}` 파일 로드 실패. GitHub에 파일이 없거나 내용이 비어있을 수 있습니다.")
            # --- ★★★ 디버깅 코드 추가 끝 ★★★ ---
        
        progress_bar.progress((i + 1) / len(smhb_codes))

    progress_bar.empty()

    if not Target_DataFrames:
        st.error("선택된 약재에 대한 유효 데이터를 불러오지 못했습니다. 위 로그를 확인해주세요.")
        return None, None, None, None

    # --- ▼▼▼ 제공해주신 분석 코드의 '의도'를 CSV(DataFrame)에 맞게 구현 ▼▼▼ ---
    
    st.info("로드된 CSV 데이터를 파싱하여 Node와 Edge 목록을 생성합니다...")
    node_list = []
    edge_list = []

    # Step 2: 각 약재의 DataFrame을 순회하며 node_list와 edge_list 생성
    # 이 부분이 CSV 형식에 맞춰 재구성된 핵심 파싱 로직입니다.
    for herb_df in Target_DataFrames:
        for i in range(len(herb_df)):
            if herb_df[i]['group'] == 'nodes' :
                node_list.append(herb_df[i]['data']['info'].split("<br>"))
            if herb_df[i]['group'] == 'edges' :
                edge_list.append([herb_df[i]['data']['source'],  herb_df[i]['data']['target']])

    # Step 3: 이하 제공해주신 코드 로직을 기반으로 데이터 가공 및 분석 수행
    # data shaping
    node_data = pd.DataFrame()
    node_ID, node_label, node_group = [], [], []
    SourceID = []
    TargetID = []
    for i in range(len(node_list)) :
        node_ID.append(node_list[i][1][4:].strip())
        node_label.append(node_list[i][2][6:])
        node_group.append(node_list[i][0])
        
    node_data['ID'] = node_ID
    node_data['Label'] =  node_label
    node_data['Group'] = node_group
    
    for i in range(len(edge_list)) :
        SourceID.append(edge_list[i][0])
        TargetID.append(edge_list[i][1])
        
    edge_data['SourceID'] = SourceID
    edge_data['TargetID'] = TargetID
    #herb origin info
    origin_list = []
    
    for i in range(len(node_data)) :
        if node_data.iloc[i]['ID'][:4] == 'SMHB' :
            name = node_data.iloc[i]['ID']
        origin_list.append(name)
    
    node_data['origin'] = origin_list
    
    #OB score loading
    ingre_data = load_excel_data('SMIT.xlsx')
    
    #OB score null data deletion 
    ingre_data = ingre_data.dropna(subset=['OB_score'], how='any', axis=0)
    
    
    
    #OB score가 없는 ingredient데이터들 제거
    drop_list = []
    
    for i in tqdm(range(len(node_data))) : 
        if node_data.iloc[i]['Group'] == 'MM symptom' : #불필요데이터제거
            drop_list.append(i)
        if node_data.iloc[i]['Group'] == 'TCM symptom' : #불필요데이터제거
            drop_list.append(i)
        if node_data.iloc[i]['Group'] != 'ingredient' : 
            continue
        if node_data.iloc[i]['Group'] == 'ingredient' :
            if node_data.iloc[i]['Label'] not in  list(ingre_data['Molecule_name']) :
                drop_list.append(i)

    st.info("OB_Score_ clear")
    # 1단계: 모든 파일에서 p-value 수집하기
    all_p_values_data = []
    print("1/3: 모든 파일에서 P-value 수집 중...")
    for p in smhb_codes:
        csv_data = load_herb_csv_data(p)
        
        # 파일 원본(origin) 정보 추가
        csv_data['origin'] = p[:9] 
        
        all_p_values_data.append(csv_data[['origin', 'Target id', 'P_value']])
    
    # 하나의 데이터프레임으로 합치기
    all_tests_df = pd.concat(all_p_values_data, ignore_index=True)
    all_tests_df.dropna(subset=['P_value'], inplace=True) # P_value가 없는 경우 제거
    
    # 2단계: FDR 보정 적용하기
    print("2/3: FDR 보정 적용 중...")
    # multipletests 함수는 p-value 리스트를 받아 q-value 리스트를 반환합니다.
    # rejected: (True/False), q_values: 보정된 p-value, ...
    rejected, q_values, _, _ = multipletests(all_tests_df['P_value'], alpha=0.05, method='fdr_bh')
    all_tests_df['q_value'] = q_values
    
    # 3단계: 유의미한 타겟 목록 만들기 (q-value < 0.05)
    significant_df = all_tests_df[all_tests_df['q_value'] < 0.05]
    
    # 빠른 조회를 위해 (origin, Target id) 쌍을 set 형태로 변환
    significant_targets_set = set(zip(significant_df['origin'], significant_df['Target id']))
    print(f"총 {len(all_tests_df)}개의 테스트 중 {len(significant_targets_set)}개의 유의미한 타겟을 발견했습니다.")
    
    
    # 4단계: 최종 노드 데이터 필터링하기
    print("3/3: 유의미하지 않은 Target 노드 제거 중...")
    # 기존 drop_list는 그대로 사용하거나, 여기서 새로 만들어도 됩니다.
    # drop_list = [] # 새로 시작할 경우 주석 해제
    
    # node_data의 target 들 중, significant_targets_set에 없는 것들을 drop_list에 추가
    for index, row in tqdm(node_data.iterrows(), total=len(node_data)):
        if row['Group'] == 'target':
            # (origin, ID) 쌍이 유의미한 목록에 없으면 제거 대상에 추가
            if (row['origin'], row['ID']) not in significant_targets_set:
                drop_list.append(index)
    
    # 중복된 인덱스 제거 후 최종적으로 노드 제거
    final_drop_indices = sorted(list(set(drop_list)))
    print(f"{len(final_drop_indices)}개의 노드를 제거합니다.")
    node_data.drop(final_drop_indices, axis=0, inplace=True)
    
    st.info(f"유효성분 필터링 완료. {len(drop_list)}개의 노드 제거.")
    
    
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
