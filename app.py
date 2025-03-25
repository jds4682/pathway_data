import streamlit as st
import networkx as nx
import pandas as pd
import os
import plotly.graph_objects as go
import requests
from io import BytesIO
from streamlit_plotly_events import plotly_events

def load_pathway_data(name):
    url = f"https://raw.githubusercontent.com/jds4682/pathway_data/main/{name}_pathway_scores.xlsx"
    response = requests.get(url)
    print("새로 url 받았음")
    if response.status_code == 200:
        return pd.read_excel(BytesIO(response.content))
        
    else:
        st.error("Could not able to download contents.")
        return None

# Define filters
T1 = ["Daeshiho-tang", "SMHB00058", "SMHB00188", "SMHB00021", "SMHB00085", "SMHB00486", "SMHB00035"]
T2 = ["Banhabakchulcheonmatang", "SMHB00035", "SMHB00064", "SMHB00272", "SMHB00022", "SMHB00052", "SMHB00336", "SMHB00187", "SMHB00527", "SMHB00129", "SMHB00480", "SMHB00136", "SMHB00184", "SMHB00367"]
T3 = ["Banhasasim-tang", "SMHB00035", "SMHB00188", "SMHB00336", "SMHB00133", "SMHB00136", "SMHB00186", "SMHB00367", "SMHB00090"]
T4 = ["Banhahubak-tang", "SMHB00035", "SMHB00129", "SMHB00171", "SMHB00497", "SMHB00367", "SMHB00090"]
T5 = ["Bojungikgi-tang", "SMHB00187", "SMHB00336", "SMHB00022", "SMHB00133", "SMHB00095", "SMHB00064", "SMHB00366", "SMHB00058"]
T6 = ["Saengmaek-san", "SMHB00336", "SMHB00041"]
T7 = ["Sosiho-tang", "SMHB00058", "SMHB00188", "SMHB00336", "SMHB00035", "SMHB00133", "SMHB00367", "SMHB00090"]
T8 = ["Socheongryong-tang", "SMHB00264", "SMHB00021", "SMHB00035", "SMHB00041", "SMHB00155", "SMHB00425", "SMHB00133", "SMHB00136"]
T9 = ["Hyeonggaeyeongyo-tang", "SMHB00212", "SMHB00241", "SMHB00122", "SMHB00095", "SMHB00077", "SMHB00058", "SMHB00485", "SMHB00188", "SMHB00027", "SMHB00216", "SMHB00021", "SMHB00483", "SMHB00133"]
T10 = ["Hwanglyeonhaedok-tang", "SMHB00186", "SMHB00188", "SMHB00184", "SMHB00483"]
T11 = ["Gamisoyo-san", "SMHB00483", "SMHB00284", "SMHB00133", "SMHB00129", "SMHB00095", "SMHB00058", "SMHB00037", "SMHB00022", "SMHB00021"]
T12 = ["Galgeun-tang", "SMHB00021", "SMHB00090", "SMHB00133", "SMHB00155", "SMHB00264", "SMHB00367", "SMHB00499"]
T13 = ["Galgeunhaegi-tang", "SMHB00021", "SMHB00027", "SMHB00058", "SMHB00090", "SMHB00133", "SMHB00188", "SMHB00216", "SMHB00322", "SMHB00366", "SMHB00367", "SMHB00499"]


T1_weights = {
    "SMHB00058": 15,
    "SMHB00188": 9.38,
    "SMHB00021": 9.38,
    "SMHB00085": 7.5,
    "SMHB00486": 5.63,
    "SMHB00035": 3.75
}

T2_weights = {
    "SMHB00035": 5.63,
    "SMHB00064": 5.63,
    "SMHB00272": 5.63,
    "SMHB00022": 3.75,
    "SMHB00052": 1.88,
    "SMHB00336": 1.88,
    "SMHB00187": 1.88,
    "SMHB00527": 1.88,
    "SMHB00129": 1.88,
    "SMHB00480": 1.88,
    "SMHB00136": 1.13,
    "SMHB00184": 0.75,
    "SMHB00367": 2.5
}

T3_weights = {
    "SMHB00035": 7.5,
    "SMHB00188": 5.63,
    "SMHB00336": 5.63,
    "SMHB00133": 5.63,
    "SMHB00136": 3.75,
    "SMHB00186": 1.88,
    "SMHB00367": 1.5,
    "SMHB00090": 2
}

T4_weights = {
    "SMHB00035": 7.5,
    "SMHB00129": 6,
    "SMHB00171": 4.5,
    "SMHB00497": 3,
    "SMHB00367": 3.5,
    "SMHB00090": 2
}

T5_weights = {
    "SMHB00187": 5.63,
    "SMHB00336": 3.75,
    "SMHB00022": 3.75,
    "SMHB00133": 3.75,
    "SMHB00095": 1.88,
    "SMHB00064": 1.88,
    "SMHB00366": 1.13,
    "SMHB00058": 1.13
}

T6_weights = {
    
    "SMHB00336" : 3.75,
    "SMHB00041" : 3.75 
}

T7_weights = {
    "SMHB00058": 11.25,
    "SMHB00188": 7.5,
    "SMHB00336": 3.75,
    "SMHB00035": 3.75,
    "SMHB00133": 1.88,
    "SMHB00367": 1.5,
    "SMHB00090": 2
}

T8_weights = {
    "SMHB00264": 5.63,
    "SMHB00021": 5.63,
    "SMHB00035": 5.63,
    "SMHB00041": 5.63,
    "SMHB00155": 3.75,
    "SMHB00425": 3.75,
    "SMHB00133": 3.75,
    "SMHB00136": 3.75
}

T9_weights = {
    "SMHB00212": 2.63,
    "SMHB00241": 2.63,
    "SMHB00122": 2.63,
    "SMHB00095": 2.63,
    "SMHB00077": 2.63,
    "SMHB00058": 2.63,
    "SMHB00485": 2.63,
    "SMHB00188": 2.63,
    "SMHB00027": 2.63,
    "SMHB00216": 2.63,
    "SMHB00021": 2.63,
    "SMHB00483": 2.63,
    "SMHB00133": 1.88
}

T10_weights = {
    "SMHB00186": 4.69,
    "SMHB00188": 4.69,
    "SMHB00184": 4.69,
    "SMHB00483": 4.69
}

T11_weights = {
    "SMHB00095": 3.75,
    "SMHB00021": 3.75,
    "SMHB00129": 3.75,
    "SMHB00022": 3.75,
    "SMHB00058": 3.75,
    "SMHB00483": 3.75,
    "SMHB00284": 3.75,
    "SMHB00133": 3.75,
    "SMHB00037": 1.88
}

T12_weights = {
    "SMHB00499": 11.25,
    "SMHB00264": 7.5,
    "SMHB00021": 5.63,
    "SMHB00155": 3.75,
    "SMHB00133": 3,
    "SMHB00367": 1.5,
    "SMHB00090": 2
}

T13_weights = {
    "SMHB00499": 3.75,
    "SMHB00058": 3.75,
    "SMHB00188": 3.75,
    "SMHB00322": 3.75,
    "SMHB00542": 3.75,
    "SMHB00021": 3.75,
    "SMHB00366": 3.75,
    "SMHB00027": 3.75,
    "SMHB00216": 3.75,
    "SMHB00133": 1.88,
    "SMHB00367": 1.5,
    "SMHB00090": 2
}
#don't touch above this area!!!!!!!!!!!!!


filters = {
    T1[0]: (T1, T1_weights),
    T2[0]: (T2, T2_weights),
    T3[0]: (T3, T3_weights),
    T4[0]: (T4, T4_weights),
    T5[0]: (T5, T5_weights),
    T6[0]: (T6, T6_weights),
    T7[0]: (T7, T7_weights),
    T8[0]: (T8, T8_weights),
    T9[0]: (T9, T9_weights),
    T10[0]: (T10, T10_weights),
    T11[0]: (T11, T11_weights),
    T12[0]: (T12, T12_weights),
    T13[0]: (T13, T13_weights),
}
filter_options = list(filters.keys())
selected_filter = st.selectbox("Select a Filter", filter_options)
selected_tang, selected_weights = filters[selected_filter]
tang_name = selected_tang[0]
st.write(tang_name)

df_pathway = load_pathway_data(tang_name)
data_list = []


G = nx.Graph()
G.add_node(tang_name, type='prescription', color='red', layer=0, size=12)

for herb in selected_tang[1:]:
    G.add_node(herb, type='herb', color='orange', layer=1, size=8)
    G.add_edge(tang_name, herb, weight=1.5)

    url_path = f"https://github.com/jds4682/pathway_data/raw/main/{herb}.csv"
    response = requests.get(url_path)

    if response.status_code == 200:
        try:
            df = pd.read_csv(BytesIO(response.content), encoding='ISO-8859-1')
        except UnicodeDecodeError:
            continue
    else:
        continue

    df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
    df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]
    filtered_df = df[(df['P_value'].astype(float) < 0.01) & (df['Value'].astype(float) > 1)]

    for _, row in filtered_df.iterrows():
        gene = row['Gene symbol']
        value = float(row['Value'])
        score = value * selected_weights.get(herb, 1)
        data_list.append([herb, gene, score])
        G.add_node(gene, type='gene', size=max(15, min(score * 0.5, 3)), color='green', layer=2)
        G.add_edge(herb, gene, weight=score)

if df_pathway is not None:
    for _, row in df_pathway.iterrows():
        gene = row['Gene']
        pathway = row['Pathway']
        score = row['Score']
        total_score = row['Total Score']

        G.add_node(pathway, type='pathway', size=max(15, min(total_score * 0.5, 3)), color='purple', layer=3)
        G.add_edge(gene, pathway, weight=score)

pathway_options = ["All"] + (list(df_pathway['Pathway'].unique()) if df_pathway is not None else [])
if 'selected_node' not in st.session_state:
    st.session_state['selected_node'] = None

pathway_filter = st.selectbox("Select a Pathway", pathway_options, index=0)

def update_graph(pathway_filter, selected_node):
    filtered_G = G.copy()

    if selected_node:
        nodes_to_keep = set([selected_node]) | set(G.neighbors(selected_node))

        # 선택한 노드의 이웃의 이웃까지 포함
        for neighbor in list(G.neighbors(selected_node)):
            nodes_to_keep.update(G.neighbors(neighbor))

        # 확장된 서브그래프 생성
        filtered_G = G.subgraph(nodes_to_keep)
    
    elif pathway_filter != "All":
        nodes_to_keep = {n for n, d in G.nodes(data=True) if d['type'] in ['prescription', 'herb']}

        for edge in G.edges():
            if edge[1] == pathway_filter or edge[0] == pathway_filter:
                nodes_to_keep.add(edge[0])
                nodes_to_keep.add(edge[1])

        filtered_G = G.subgraph(nodes_to_keep)

    # 필터링된 그래프의 레이아웃 설정
    layers = {'prescription': [], 'herb': [], 'gene': [], 'pathway': []}
    for node, data in filtered_G.nodes(data=True):
        layers[data['type']].append(node)

    shell_positions = nx.shell_layout(
        filtered_G, 
        [layers['prescription'], layers['herb'], layers['gene'], layers['pathway']]
    )

    nodes = list(filtered_G.nodes())
    node_colors = [filtered_G.nodes[n]['color'] for n in nodes]
    node_sizes = [filtered_G.nodes[n].get('size', 300) for n in nodes]
    edges = list(filtered_G.edges())

    edge_x, edge_y = [], []
    for edge in edges:
        x0, y0 = shell_positions[edge[0]]
        x1, y1 = shell_positions[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1, color='gray'),
        hoverinfo='none',
        mode='lines',
        name='Edges'
    ))

    for node_type, color in [('prescription', 'red'), ('herb', 'orange'), ('gene', 'green'), ('pathway', 'purple')]:
        filtered_nodes = [n for n in nodes if filtered_G.nodes[n]['color'] == color]
        if filtered_nodes:
            fig.add_trace(go.Scatter(
                x=[shell_positions[n][0] for n in filtered_nodes],
                y=[shell_positions[n][1] for n in filtered_nodes],
                mode='markers+text',
                text=filtered_nodes,
                marker=dict(size=[filtered_G.nodes[n].get('size', 300) for n in filtered_nodes], color=color, opacity=0.8),
                textposition='top center',
                name=node_type
            ))

    fig.update_layout(
        width=1200,
        height=1000,
        showlegend=True,
        title=f"{tang_name} Network Graph",
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False)
    )

    fig.update_traces(marker=dict(symbol='circle'), selector=dict(mode='markers+text'))
    fig.update_layout(clickmode='event+select')

    # 선택된 노드 이벤트 가져오기
    
    
    return fig
fig = update_graph(pathway_filter, st.session_state['selected_node'])
# ✅ 클릭 이벤트 추가
selected_points = plotly_events(fig)  # 사용자가 클릭한 포인트 가져오기

# ✅ 선택된 노드 업데이트
if selected_points:
    point_index = selected_points[0]["pointIndex"]  # 클릭된 점의 index 가져오기
    nodes = list(G.nodes())  # 현재 그래프의 노드 리스트 가져오기
    
    if point_index < len(nodes):  # 인덱스가 범위를 벗어나지 않는지 확인
        new_selected_node = nodes[point_index]  # 새로운 선택된 노드
        
        if st.session_state.get("selected_node") != new_selected_node:  # 기존 노드와 다를 때만 갱신
            st.session_state["selected_node"] = new_selected_node  # 선택된 노드 저장
            st.rerun()  # ✅ 새 그래프를 갱신하기 위해 페이지 다시 로드

# ✅ 선택된 노드 표시
if "selected_node" in st.session_state and st.session_state["selected_node"]:
    st.write(f"선택된 노드: {st.session_state['selected_node']}")

# ✅ 선택된 노드가 반영된 그래프 생성
fig = update_graph(pathway_filter, st.session_state.get("selected_node"))
st.plotly_chart(fig)

# Reset 버튼 추가
if st.button("Reset Selection"):
    st.session_state["selected_node"] = None
    st.rerun()
