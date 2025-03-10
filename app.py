import streamlit as st
import networkx as nx
import pandas as pd
import os
import plotly.graph_objects as go
import requests
from io import BytesIO

url = "https://github.com/jds4682/pathway_data/raw/db346c9671fd44fc808ffe11cbc3b2bc788513d9/Saengmaek-san_pathway_scores.xlsx"
response = requests.get(url)

if response.status_code == 200:
    df_pathway = pd.read_excel(BytesIO(response.content))
else:
    st.error("파일을 다운로드할 수 없습니다.")

T6 = ["Saengmaek-san", "SMHB00336", "SMHB00041"]
T6_weights = {"SMHB00336": 3.75, "SMHB00041": 3.75}

data_list = []
tang_name = T6[0]
st.write(tang_name)

G = nx.Graph()
G.add_node(tang_name, type='prescription', color='red', layer=0, size=12)

for herb in T6[1:]:
    G.add_node(herb, type='herb', color='orange', layer=1, size=8)
    G.add_edge(tang_name, herb, weight=1.5)
    
    url_path = f"https://github.com/jds4682/pathway_data/raw/refs/heads/main/{herb}.csv"
    response = requests.get(url_path)
    
    if response.status_code == 200:
        df = pd.read_csv(BytesIO(response.content), encoding='ISO-8859-1')
    else:
        continue
    
    df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
    df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]
    
    filtered_df = df[(df['P_value'].astype(float) < 0.01) & (df['Value'].astype(float) > 1)]
    
    for _, row in filtered_df.iterrows():
        gene = row['Gene symbol']
        value = float(row['Value'])
        score = value * T6_weights.get(herb, 1)
        G.add_node(gene, type='gene', size=max(15, min(score * 0.5, 3)), color='green', layer=2)
        G.add_edge(herb, gene, weight=score)

pathway_options = ["All"] + list(df_pathway['Pathway'].unique())
pathway_filter = st.selectbox("Select a Pathway", pathway_options, index=0)

selected_node = st.session_state.get("selected_node", None)

def update_graph(pathway_filter, selected_node):
    filtered_G = G.copy()
    if selected_node:
        nodes_to_keep = {selected_node} | {neighbor for neighbor in G.neighbors(selected_node)}
        filtered_G = G.subgraph(nodes_to_keep)
    elif pathway_filter != "All":
        nodes_to_keep = {n for n, d in G.nodes(data=True) if d['type'] in ['prescription', 'herb']}
        for edge in G.edges():
            if edge[1] == pathway_filter or edge[0] == pathway_filter:
                nodes_to_keep.add(edge[0])
                nodes_to_keep.add(edge[1])
        filtered_G = G.subgraph(nodes_to_keep)
    
    pos = nx.spring_layout(G, seed=42)
    nodes = list(filtered_G.nodes())
    node_colors = [filtered_G.nodes[n]['color'] for n in nodes]
    node_sizes = [filtered_G.nodes[n].get('size', 300) for n in nodes]
    edges = list(filtered_G.edges())
    edge_x, edge_y = [], []
    for edge in edges:
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1, color='gray'),
        hoverinfo='none',
        mode='lines'))
    
    fig.add_trace(go.Scatter(
        x=[pos[n][0] for n in nodes],
        y=[pos[n][1] for n in nodes],
        mode='markers+text',
        text=nodes,
        marker=dict(size=node_sizes, color=node_colors, opacity=0.8),
        textposition='top center',
        customdata=nodes,
        hoverinfo='text'))
    
    fig.update_layout(
        width=1200,
        height=1000,
        showlegend=False,
        title=f"{tang_name} Network Graph",
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        clickmode='event+select')
    
    return fig

fig = update_graph(pathway_filter, selected_node)
st.plotly_chart(fig)

if st.session_state.get("selected_node", None):
    st.write(f"현재 선택된 노드: {st.session_state.selected_node}")
