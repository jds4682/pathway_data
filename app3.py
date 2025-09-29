# 필요한 추가 import (맨 위에 이미 있으면 중복 안 해도 됩니다)
import gseapy as gp
import mygene
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os
from shutil import copyfile

def process_and_run_gsea_gseapy_fix(prescription_name, selected_herbs_info, herb_weights):
    # 1) 기존 전처리 (same as before) -> py_df
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

    # 2) Aggregate by symbol
    agg = py_df.groupby('GeneSymbol', as_index=False).agg({'Score': 'sum'})
    agg['GeneSymbol'] = agg['GeneSymbol'].astype(str).str.strip()

    # 3) mygene 매핑: symbol 우선, alias도 탐색
    mg = mygene.MyGeneInfo()
    symbols = agg['GeneSymbol'].tolist()
    try:
        mg_res = mg.querymany(symbols, scopes='symbol,alias', fields='symbol,entrezgene', species='human', as_dataframe=True)
    except Exception as e:
        st.warning(f"mygene querymany 에러: {e}. SYMBOL 기반으로 그대로 진행합니다.")
        mg_res = None

    if isinstance(mg_res, pd.DataFrame):
        # mg_res index is the original query
        mg_res = mg_res.reset_index().rename(columns={'query': 'query_symbol'}) if 'query' in mg_res.columns else mg_res.reset_index().rename(columns={'index':'query_symbol'})
        # prefer canonical 'symbol' returned by mygene; if missing, keep original
        mapping_df = mg_res[['query_symbol', 'symbol']].drop_duplicates(subset=['query_symbol'])
        merged = pd.merge(agg, mapping_df, left_on='GeneSymbol', right_on='query_symbol', how='left')
        merged['gsea_id'] = merged['symbol'].fillna(merged['GeneSymbol']).astype(str).str.upper()
        # unmapped
        unmapped = merged[merged['symbol'].isnull()]['GeneSymbol'].unique().tolist()
    else:
        merged = agg.copy()
        merged['gsea_id'] = merged['GeneSymbol'].astype(str).str.upper()
        unmapped = []

    st.info(f"Input genes: {len(agg)}, unmapped (example up to 20): {unmapped[:20]}")
    # 4) Build ranked Series (symbol uppercase)
    ranked = merged.groupby('gsea_id', as_index=False).agg({'Score':'sum'}).sort_values('Score', ascending=False)
    ranked = ranked.set_index('gsea_id')['Score'].astype(float)

    # 5) Tie-breaking: deterministic tiny jitter to make all values unique
    eps = 1e-9
    ranked = ranked + (np.arange(len(ranked)) * eps)

    # 6) (디버깅) gene-set overlap 확인 (예: GO_Biological_Process_2021)
    try:
        # download library dict (internal gseapy helper)
        from gseapy.parser import download_library
        lib_name = 'GO_Biological_Process_2021'
        libdict = download_library(lib_name, organism='Human')
        lib_genes = set(sum(libdict.values(), []))
        overlap = set(ranked.index) & lib_genes
        st.info(f"Overlap with {lib_name}: {len(overlap)} genes (입력 총 {len(ranked)})")
    except Exception as e:
        st.info(f"라이브러리 overlap 체크 실패 (네트워크/내부함수 문제): {e}")

    # 7) Run prerank with threads and relaxed min_size (필요시 조정)
    tmpdir = tempfile.mkdtemp(prefix='gsea_out_')
    out_plots = {}
    try:
        pre_res_go = gp.prerank(rnk=ranked, gene_sets='GO_Biological_Process_2021',
                                threads=4, permutation_num=500, outdir=tmpdir,
                                format='png', seed=123, min_size=5, max_size=2000)
        res_go = pre_res_go.res2d
        if res_go is not None and not res_go.empty:
            # 간단한 dotplot (예시)
            df_top = res_go.reset_index().head(15)
            df_top['neglog10fdr'] = -np.log10(df_top['fdr_q-val'].replace(0, np.nextafter(0,1)))
            plt.figure(figsize=(10,6))
            sizes = df_top.get('geneset_size', np.ones(len(df_top))).astype(float)
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
            # gseaplot for top (if available)
            try:
                top_term = res_go.index[0]
                gp.plot.gseaplot(pre_res_go.ranking, pre_res_go.results[top_term], ofname=os.path.join(tmpdir, 'plot_go_gseaplot.png'))
                out_plots['plot_go_gseaplot.png'] = os.path.join(tmpdir, 'plot_go_gseaplot.png')
            except Exception:
                pass
    except Exception as e:
        st.warning(f"GO prerank 수행 중 오류: {e}. (min_size 조정 필요할 수 있음)")

    # 동일 방식으로 KEGG 실행 (필요 시)
    try:
        pre_res_kegg = gp.prerank(rnk=ranked, gene_sets='KEGG_2019_Human',
                                  threads=4, permutation_num=500, outdir=tmpdir,
                                  format='png', seed=123, min_size=5, max_size=2000)
        res_kegg = pre_res_kegg.res2d
        if res_kegg is not None and not res_kegg.empty:
            df_top = res_kegg.reset_index().head(15)
            df_top['neglog10fdr'] = -np.log10(df_top['fdr_q-val'].replace(0, np.nextafter(0,1)))
            plt.figure(figsize=(10,6))
            sizes = df_top.get('geneset_size', np.ones(len(df_top))).astype(float)
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
            try:
                top_term = res_kegg.index[0]
                gp.plot.gseaplot(pre_res_kegg.ranking, pre_res_kegg.results[top_term], ofname=os.path.join(tmpdir, 'plot_kegg_gseaplot.png'))
                out_plots['plot_kegg_gseaplot.png'] = os.path.join(tmpdir, 'plot_kegg_gseaplot.png')
            except Exception:
                pass
    except Exception as e:
        st.warning(f"KEGG prerank 수행 중 오류: {e}. (min_size 조정 필요할 수 있음)")

    if out_plots:
        # 복사하여 반환
        final = {}
        persistent_dir = tempfile.mkdtemp(prefix='gsea_final_')
        for name, p in out_plots.items():
            if os.path.exists(p):
                dest = os.path.join(persistent_dir, name)
                copyfile(p, dest)
                final[name] = dest
        st.success("GSEA (python) 완료")
        return final
    else:
        st.error("분석은 완료되었으나 저장된 플롯이 없습니다. (매핑/파라미터 확인 필요)")
        return None
