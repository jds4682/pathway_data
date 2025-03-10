import streamlit as st
import networkx as nx
import pandas as pd
import os
import plotly.graph_objects as go

import requests
import pandas as pd
from io import BytesIO

url = "https://github.com/jds4682/pathway_data/raw/db346c9671fd44fc808ffe11cbc3b2bc788513d9/Saengmaek-san_pathway_scores.xlsx"
response = requests.get(url)

# 파일 내용 읽기
if response.status_code == 200:
    df_pathway = pd.read_excel(BytesIO(response.content))
    print(df_pathway.head())
else:
    print("파일을 다운로드할 수 없습니다.")


# 데이터 파일 경로 설정
path = r'C:\tg'  # 데이터 파일이 있는 경로

T6 = ["Saengmaek-san", "SMHB00336", "SMHB00041"]
T6_weights = {
    "SMHB00336": 3.75,
    "SMHB00041": 3.75
}

data_list = []

tang_name = T6[0]
st.write(tang_name)

# 네트워크 그래프 생성
G = nx.Graph()
G.add_node(tang_name, type='prescription', color='red', layer=0, size=12)

for herb in T6[1:]:
    G.add_node(herb, type='herb', color='orange', layer=1, size=8)
    G.add_edge(tang_name, herb, weight=1.5)
    
    
    url_path = "https://github.com/jds4682/pathway_data/raw/refs/heads/main/" + herb + '.csv'
    print(url_path)
    response = requests.get(url_path)
    
    # 파일 내용 읽기
    if response.status_code == 200:
        try:
            df = pd.read_csv(BytesIO(response.content), encoding='ISO-8859-1')  # 다른 인코딩 시도
            print(df.head())
        except UnicodeDecodeError:
            print("파일을 읽을 수 없습니다. 인코딩 문제 발생.")
    else:
        print("파일을 다운로드할 수 없습니다.")
    
    
    
    df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
    df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]
    
    filtered_df = df[(df['P_value'].astype(float) < 0.01) & (df['Value'].astype(float) > 1)]
    
    for _, row in filtered_df.iterrows():
        gene = row['Gene symbol']
        value = float(row['Value'])
        score = value * T6_weights.get(herb, 1)
        data_list.append([herb, gene, score])
        G.add_node(gene, type='gene', size=max(15, min(score * 0.5, 3)), color='green', layer=2)
        G.add_edge(herb, gene, weight=score)

# Pathway 데이터 불러오기
pathway_file = 1
if os.path.exists(pathway_file):
    
    for _, row in df_pathway.iterrows():
        gene = row['Gene']
        pathway = row['Pathway']
        score = row['Score']
        total_score = row['Total Score']
        
        G.add_node(pathway, type='pathway', size=max(15, min(total_score * 0.5, 3)), color='purple', layer=3)
        G.add_edge(gene, pathway, weight=score)

pathway_options = ["All"] + list(df_pathway['Pathway'].unique())

# Streamlit 애플리케이션 레이아웃
pathway_filter = st.selectbox("Select a Pathway", pathway_options, index=0)

# 그래프 업데이트
def update_graph(pathway_filter):
    filtered_G = G.copy()
    if pathway_filter != "All":
        nodes_to_keep = {n for n, d in G.nodes(data=True) if d['type'] in ['prescription', 'herb']}
        for edge in G.edges():
            if edge[1] == pathway_filter or edge[0] == pathway_filter:
                nodes_to_keep.add(edge[0])
                nodes_to_keep.add(edge[1])
        filtered_G = G.subgraph(nodes_to_keep)
    
    layers = {'prescription': [], 'herb': [], 'gene': [], 'pathway': []}
    for node, data in filtered_G.nodes(data=True):
        layers[data['type']].append(node)
    
    shell_positions = nx.shell_layout(filtered_G, 
        [layers['prescription'], layers['herb'], layers['gene'], layers['pathway']])
    
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
        name='Edges'))
    
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
                name=node_type))
    
    fig.update_layout(
        width=1200,
        height=1000,
        showlegend=True,
        title=f"{tang_name} Network Graph",
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False))
    
    # HTML 파일로 내보내기
    fig.write_html(f"{tang_name}_network_graph.html")
    
    return fig

# 선택된 Pathway에 따라 그래프 업데이트
fig = update_graph(pathway_filter)

# Streamlit에 그래프 표시
st.plotly_chart(fig)
