import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import os
import tempfile

import gseapy as gp
import mygene
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tempfile
import os

# --- 대체할 함수: R 대신 gseapy 사용 ---
def process_and_run_gsea_gseapy(prescription_name, selected_herbs_info, herb_weights):
    # Step 1: 기존과 동일한 전처리 -> py_df 생성
    data_list = []
    for herb_name, herb_code in selected_herbs_info.items():
        df = load_herb_csv_data(herb_code)
        if df is None or df.empty:
            continue
        df = df[pd.to_numeric(df['P_value'], errors='coerce').notna()]
        df = df[pd.to_numeric(df['Value'], errors='coerce').notna()]
        weight = herb_weights.get(herb_name, 1.0)
        for _, row in df.iterrows():
            data_list.append([herb_code, row['Gene symbol'], float(row['Value']) * weight])

    if not data_list:
        st.error("처리할 데이터가 없습니다.")
        return None

    py_df = pd.DataFrame(data_list, columns=['herb', 'GeneSymbol', 'Score'])

    # aggregate by GeneSymbol (like clusterProfiler step)
    agg = py_df.groupby('GeneSymbol', as_index=False).agg({'Score': 'sum'})
    agg = agg.dropna(subset=['GeneSymbol'])
    if agg.empty:
        st.error("유효한 유전자 데이터가 없습니다.")
        return None

    # Step 2: SYMBOL -> ENTREZID 매핑 (mygene 사용)
    mg = mygene.MyGeneInfo()
    symbols = agg['GeneSymbol'].astype(str).tolist()
    try:
        mapping = mg.querymany(symbols, scopes='symbol', fields='entrezgene', species='human', as_dataframe=True)
    except Exception as e:
        st.warning(f"Gene ID 매핑 중 오류 발생 (mygene): {e}. SYMBOL 기반으로 진행합니다.")
        mapping = None

    if isinstance(mapping, pd.DataFrame) and 'entrezgene' in mapping.columns:
        # querymany returns index = query symbol; ensure alignment
        mapping = mapping.reset_index().rename(columns={'index': 'query'})
        mapping = mapping[['query', 'entrezgene']].drop_duplicates(subset=['query'])
        merged = pd.merge(agg, mapping, left_on='GeneSymbol', right_on='query', how='left')
        # prefer entrez for ranking if available, else keep SYMBOL
        merged['gene_id_for_gsea'] = merged['entrezgene'].fillna(merged['GeneSymbol'])
    else:
        merged = agg.copy()
        merged['gene_id_for_gsea'] = merged['GeneSymbol']

    # create preranked list: Series index = gene name (entrez or symbol), value = Score
    # gseapy expects a file or pd.Series with gene_name and rank
    ranked = merged[['gene_id_for_gsea', 'Score']].copy()
    ranked = ranked.dropna(subset=['gene_id_for_gsea'])
    # convert entrez floats to str, keep unique by taking sum if duplicates
    ranked['gene_id_for_gsea'] = ranked['gene_id_for_gsea'].astype(str)
    ranked = ranked.groupby('gene_id_for_gsea', as_index=False).agg({'Score': 'sum'})
    ranked = ranked.sort_values(by='Score', ascending=False)

    # prepare temp dir to save plots
    with tempfile.TemporaryDirectory() as tmpdir:
        out_plots = {}

        # 1) GSEA for GO (using Enrichr library name for GO BP)
        try:
            pre_res_go = gp.prerank(rnk=ranked, gene_sets='GO_Biological_Process_2021',
                                    processes=4, permutation_num=1000, outdir=tmpdir,
                                    format='png', seed=123)
            # pre_res_go.res2d is the results DataFrame
            res_go = pre_res_go.res2d
            if res_go is not None and not res_go.empty:
                # save a dotplot-like figure: NES vs -log10(fdr)
                df_top = res_go.reset_index().head(15)
                df_top['neglog10fdr'] = -np.log10(df_top['fdr_q-val'].replace(0, np.nextafter(0,1)))
                plt.figure(figsize=(10,6))
                sizes = (df_top['geneset_size'] if 'geneset_size' in df_top.columns else df_top['geneset_size']).astype(float)
                plt.scatter(df_top['NES'], df_top['neglog10fdr'], s=(sizes/np.max(sizes))*300 + 20)
                for i, txt in enumerate(df_top.index):
                    plt.text(df_top['NES'].iloc[i], df_top['neglog10fdr'].iloc[i], df_top['Term'].iloc[i], fontsize=8)
                plt.xlabel('NES')
                plt.ylabel('-log10(FDR)')
                plt.title('GO (BP) top terms')
                go_dot = os.path.join(tmpdir, 'plot_go_dotplot.png')
                plt.tight_layout()
                plt.savefig(go_dot, dpi=150)
                plt.close()
                out_plots['plot_go_dotplot.png'] = go_dot

                # GSEA plot for top 1 (if available) using gseapy's gseaplot
                try:
                    top_term = res_go.index[0]
                    gp.plot.gseaplot(pre_res_go.ranking, pre_res_go.results[top_term], ofname=os.path.join(tmpdir, 'plot_go_gseaplot.png'))
                    out_plots['plot_go_gseaplot.png'] = os.path.join(tmpdir, 'plot_go_gseaplot.png')
                except Exception:
                    pass
        except Exception as e:
            st.warning(f"GO GSEA 수행 중 오류: {e}")

        # 2) GSEA for KEGG
        try:
            pre_res_kegg = gp.prerank(rnk=ranked, gene_sets='KEGG_2019_Human',
                                       processes=4, permutation_num=1000, outdir=tmpdir,
                                       format='png', seed=123)
            res_kegg = pre_res_kegg.res2d
            if res_kegg is not None and not res_kegg.empty:
                df_top = res_kegg.reset_index().head(15)
                df_top['neglog10fdr'] = -np.log10(df_top['fdr_q-val'].replace(0, np.nextafter(0,1)))
                plt.figure(figsize=(10,6))
                sizes = (df_top['geneset_size'] if 'geneset_size' in df_top.columns else df_top['geneset_size']).astype(float)
                plt.scatter(df_top['NES'], df_top['neglog10fdr'], s=(sizes/np.max(sizes))*300 + 20)
                for i, txt in enumerate(df_top.index):
                    plt.text(df_top['NES'].iloc[i], df_top['neglog10fdr'].iloc[i], df_top['Term'].iloc[i], fontsize=8)
                plt.xlabel('NES')
                plt.ylabel('-log10(FDR)')
                plt.title('KEGG top terms')
                kegg_dot = os.path.join(tmpdir, 'plot_kegg_dotplot.png')
                plt.tight_layout()
                plt.savefig(kegg_dot, dpi=150)
                plt.close()
                out_plots['plot_kegg_dotplot.png'] = kegg_dot

                # GSEA plot for top 1 (if available)
                try:
                    top_term = res_kegg.index[0]
                    gp.plot.gseaplot(pre_res_kegg.ranking, pre_res_kegg.results[top_term], ofname=os.path.join(tmpdir, 'plot_kegg_gseaplot.png'))
                    out_plots['plot_kegg_gseaplot.png'] = os.path.join(tmpdir, 'plot_kegg_gseaplot.png')
                except Exception:
                    pass
        except Exception as e:
            st.warning(f"KEGG GSEA 수행 중 오류: {e}")

        # 3) Return paths (copy to a persistent place if needed)
        # Since tmpdir will be removed on exit, copy to a new temp folder under /tmp and return paths
        persistent_dir = tempfile.mkdtemp(prefix='gsea_out_')
        final_plots = {}
        for name, path in out_plots.items():
            if os.path.exists(path):
                dest = os.path.join(persistent_dir, name)
                try:
                    from shutil import copyfile
                    copyfile(path, dest)
                    final_plots[name] = dest
                except Exception:
                    pass

        if final_plots:
            st.success("GSEA 분석(파이썬 기반)이 완료되었습니다.")
            return final_plots
        else:
            st.error("분석은 완료되었으나 저장된 플롯이 없습니다.")
            return None

# --- 1. 초기 설정 및 GitHub 데이터 로딩 함수 ---

st.set_page_config(layout="wide", page_title="GSEA Pre-processing & Analysis Tool")

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
    return herb_df

# --- 2. GSEA 전처리 및 R 코드 실행 로직 (rpy2 수정) ---



# --- 3. 웹페이지 UI 구성 ---
st.title("🌿 GSEA 분석 자동화 웹 앱 (Docker & rpy2)")
st.info("약재와 용량을 선택하면 Python 내에서 R 코드를 직접 실행하여 GSEA 분석을 수행합니다.")

herb_df = load_initial_data()

if herb_df is not None:
    st.header("1. 처방 구성하기")
    
    KOREAN_NAME_COLUMN = 'korean name'
    SMHB_ID_COLUMN = 'SMHB_ID'
    
    try:
        herb_names = herb_df[KOREAN_NAME_COLUMN].dropna().unique().tolist()
        selected_herb_names = st.multiselect("분석에 포함할 약재를 선택하세요.", options=herb_names)
        
        selected_herbs_info = {}
        herb_weights = {}

        if selected_herb_names:
            st.subheader("용량 입력 (단위: g)")
            num_columns = 3
            cols = st.columns(num_columns)
            
            for i, name in enumerate(selected_herb_names):
                with cols[i % num_columns]:
                    grams = st.number_input(f"{name}", min_value=0.1, value=4.0, step=0.1, key=name)
                    smhb_id = herb_df[herb_df[KOREAN_NAME_COLUMN] == name][SMHB_ID_COLUMN].iloc[0]
                    selected_herbs_info[name] = smhb_id
                    herb_weights[name] = grams
            
            st.divider()
            st.header("2. 분석 실행")
            
            prescription_name_input = st.text_input("분석할 처방의 영문 이름을 입력하세요 (예: My_Prescription):")
            
            if st.button("GSEA 분석 시작", disabled=(not prescription_name_input)):
                
                plots = process_and_run_gsea_gseapy(prescription_name_input, selected_herbs_info, herb_weights)
                
                if plots:
                    st.header(f"📈 '{prescription_name_input}' GSEA 분석 결과")
                    
                    st.subheader("Gene Ontology (GO) 분석 결과")
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        if "plot_go_dotplot.png" in plots: st.image(plots["plot_go_dotplot.png"], caption="GO Dot Plot")
                    with col2:
                        if "plot_go_ridgeplot.png" in plots: st.image(plots["plot_go_ridgeplot.png"], caption="GO Ridge Plot")
                    with col3:
                        if "plot_go_gseaplot.png" in plots: st.image(plots["plot_go_gseaplot.png"], caption="GO GSEA Plot")

                    st.subheader("KEGG Pathway 분석 결과")
                    col4, col5, col6 = st.columns(3)
                    with col4:
                        if "plot_kegg_dotplot.png" in plots: st.image(plots["plot_kegg_dotplot.png"], caption="KEGG Dot Plot")
                    with col5:
                        if "plot_kegg_ridgeplot.png" in plots: st.image(plots["plot_kegg_ridgeplot.png"], caption="KEGG Ridge Plot")
                    with col6:
                        if "plot_kegg_gseaplot.png" in plots: st.image(plots["plot_kegg_gseaplot.png"], caption="KEGG GSEA Plot")

    except KeyError as e:
        st.error(f"'{e}' 열을 'all name.xlsx' 파일에서 찾을 수 없습니다. 코드의 열 이름을 확인하세요.")

