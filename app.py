import streamlit as st
import networkx as nx
import pandas as pd
import os
import plotly.graph_objects as go
import requests
from io import BytesIO

def load_pathway_data(name):
    url = "https://github.com/jds4682/pathway_data/raw/db346c9671fd44fc808ffe11cbc3b2bc788513d9/" name + "_pathway_scores.xlsx"
    response = requests.get(url)
    if response.status_code == 200:
        return pd.read_excel(BytesIO(response.content))
    else:
        st.error("파일을 다운로드할 수 없습니다.")
        return None

# Define filters
T6 = ["Saengmaek-san", "SMHB00336", "SMHB00041"]
T6_weights = {"SMHB00336": 3.75, "SMHB00041": 3.75}

T7 = ["Sosiho-tang", "SMHB00058", "SMHB00188", "SMHB00336", "SMHB00035", "SMHB00133", "SMHB00367", "SMHB00090"]
T7_weights = {
    "SMHB00058": 11.25,
    "SMHB00188": 7.5,
    "SMHB00336": 3.75,
    "SMHB00035": 3.75,
    "SMHB00133": 1.88,
    "SMHB00367": 1.5,
    "SMHB00090": 2
}

df_pathway = load_pathway_data(T6[0])
filters = {T6[0]: (T6, T6_weights), T7[0]: (T7, T7_weights)}
filter_options = list(filters.keys())
selected_filter = st.selectbox("Select a Filter", filter_options)

selected_tang, selected_weights = filters[selected_filter]

data_list = []

tang_name = selected_tang[0]
st.write(tang_name)

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
        for neighbor in list(G.neighbors(selected_node)):
            nodes_to_keep.update(G.neighbors(neighbor))
        filtered_G = G.subgraph(nodes_to_keep)
    elif pathway_filter != "All":
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

    fig.update_traces(marker=dict(symbol='circle'), selector=dict(mode='markers+text'))
    fig.update_layout(clickmode='event+select')

    return fig

fig = update_graph(pathway_filter, st.session_state['selected_node'])
st.plotly_chart(fig)  

if st.button("Reset Selection"):
    st.session_state['selected_node'] = None
    st.rerun()
