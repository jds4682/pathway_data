import streamlit as st
import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import requests
from io import BytesIO

# 데이터 불러오기 (Excel 파일에서 Pathway 정보)
url = "https://github.com/jds4682/pathway_data/raw/db346c9671fd44fc808ffe11cbc3b2bc788513d9/Saengmaek-san_pathway_scores.xlsx"
response = requests.get(url)
if response.status_code == 200:
    df_pathway = pd.read_excel(BytesIO(response.content))
else:
    st.error("Pathway 데이터를 불러올 수 없습니다.")
    df_pathway = pd.DataFrame()

# 네트워크 그래프 생성
G = nx.Graph()
T6 = ["Saengmaek-san", "SMHB00336", "SMHB00041"]
T6_weights = {"SMHB00336": 3.75, "SMHB00041": 3.75}

tang_name = T6[0]
G.add_node(tang_name, type='prescription', color='red', size=12)

data_list = []
for herb in T6[1:]:
    G.add_node(herb, type='herb', color='orange', size=8)
    G.add_edge(tang_name, herb, weight=1.5)
    
    # 유전자 데이터 불러오기
    url_path = f"https://github.com/jds4682/pathway_data/raw/refs/heads/main/{herb}.csv"
    response = requests.get(url_path)
    if response.status_code == 200:
        df = pd.read_csv(BytesIO(response.content), encoding='ISO-8859-1')
        df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
        df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]
        filtered_df = df[(df['P_value'].astype(float) < 0.01) & (df['Value'].astype(float) > 1)]
        for _, row in filtered_df.iterrows():
            gene = row['Gene symbol']
            value = float(row['Value']) * T6_weights.get(herb, 1)
            G.add_node(gene, type='gene', color='green', size=max(15, min(value * 0.5, 3)))
            G.add_edge(herb, gene, weight=value)

# Pathway 데이터 추가
for _, row in df_pathway.iterrows():
    gene, pathway, score, total_score = row['Gene'], row['Pathway'], row['Score'], row['Total Score']
    G.add_node(pathway, type='pathway', color='purple', size=max(15, min(total_score * 0.5, 3)))
    G.add_edge(gene, pathway, weight=score)

# Streamlit UI
title = st.title("Interactive Network Graph")
selected_node = st.selectbox("노드 선택", options=list(G.nodes), index=0)

# 선택한 노드와 연결된 노드만 표시하는 기능
def filter_graph(selected_node):
    nodes_to_keep = {selected_node}
    for u, v in G.edges(selected_node):
        nodes_to_keep.add(v)
    return G.subgraph(nodes_to_keep)

# 그래프 업데이트
filtered_G = filter_graph(selected_node)
def draw_graph(G):
    pos = nx.spring_layout(G)
    node_trace = go.Scatter(
        x=[pos[n][0] for n in G.nodes()],
        y=[pos[n][1] for n in G.nodes()],
        text=list(G.nodes()),
        mode='markers+text',
        marker=dict(size=[G.nodes[n].get('size', 10) for n in G.nodes()],
                    color=[G.nodes[n]['color'] for n in G.nodes()], opacity=0.8),
        textposition='top center'
    )
    edge_trace = go.Scatter(
        x=sum([[pos[u][0], pos[v][0], None] for u, v in G.edges()], []),
        y=sum([[pos[u][1], pos[v][1], None] for u, v in G.edges()], []),
        line=dict(width=1, color='gray'),
        mode='lines'
    )
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(showlegend=False, title="Filtered Network Graph",
                      xaxis=dict(showgrid=False, zeroline=False, visible=False),
                      yaxis=dict(showgrid=False, zeroline=False, visible=False))
    return fig

st.plotly_chart(draw_graph(filtered_G))
