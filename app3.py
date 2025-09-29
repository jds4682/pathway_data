import streamlit as st
import pandas as pd
import requests
from io import BytesIO
import os
import tempfile

# rpy2 관련 라이브러리 임포트
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import conversion

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

def process_and_run_gsea_rpy2(prescription_name, selected_herbs_info, herb_weights):
    
    # Step 1: Python으로 GSEA 전처리 파일 생성
    st.info("Python으로 GSEA 전처리 데이터를 생성합니다...")
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
    
    # Step 2: R 코드 실행 준비
    st.info("R 분석 환경을 설정하고 데이터를 전달합니다...")
    try:
        with conversion.localconverter(robjects.default_converter + pandas2ri.converter):
            r_df = robjects.conversion.py2rpy(py_df)

        # --- ★★★ R 코드 수정: 경로 설정 및 라이브러리 로드만 수행 ★★★ ---
        r_code = """
        # 1. Dockerfile에서 설정한 환경 변수로부터 라이브러리 경로를 가져옵니다.
        lib_path <- Sys.getenv("R_LIBS_USER")
        # 2. 이 경로를 라이브러리 검색 경로에 추가합니다.
        .libPaths(c(lib_path, .libPaths()))

        # 3. 라이브러리를 로드합니다. (이제 패키지를 찾을 수 있습니다)
        suppressPackageStartupMessages(library(clusterProfiler))
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        suppressPackageStartupMessages(library(enrichplot))
        suppressPackageStartupMessages(library(dplyr))
        suppressPackageStartupMessages(library(ggplot2))

        # 4. 메인 분석 함수 (이전과 동일)
        run_gsea_in_r <- function(gene_data_df, output_dir) {
            
            aggregated_gene_data <- gene_data_df %>%
              group_by(GeneSymbol) %>%
              summarise(TotalScore = sum(Score, na.rm = TRUE)) %>%
              as.data.frame()

            tryCatch({
                ids <- bitr(aggregated_gene_data$GeneSymbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE)
                
                gene_data_merged <- merge(aggregated_gene_data, ids, by.x="GeneSymbol", by.y="SYMBOL", all.x = TRUE)
                gene_data_final <- gene_data_merged %>% filter(!is.na(ENTREZID))
                
                geneList <- gene_data_final$TotalScore
                names(geneList) <- gene_data_final$ENTREZID
                geneList <- sort(geneList, decreasing = TRUE)
                geneList <- geneList[!duplicated(names(geneList))]
                
                if (length(geneList) == 0) {
                    print("No valid genes left after ID conversion.")
                    return()
                }

                # Run GSEA for GO
                gse_go_results <- gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont="BP", minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE, scoreType="pos")
                if (!is.null(gse_go_results) && nrow(as.data.frame(gse_go_results)) > 0) {
                    p1 <- dotplot(gse_go_results, showCategory=15)
                    ggsave(file.path(output_dir, "plot_go_dotplot.png"), plot = p1, width=10, height=8)
                    p2 <- ridgeplot(gse_go_results, showCategory=15)
                    ggsave(file.path(output_dir, "plot_go_ridgeplot.png"), plot = p2, width=10, height=8)
                    p3 <- gseaplot2(gse_go_results, geneSetID = 1:min(3, nrow(as.data.frame(gse_go_results))))
                    ggsave(file.path(output_dir, "plot_go_gseaplot.png"), plot = p3, width=10, height=8)
                }

                # Run GSEA for KEGG
                gse_kegg_results <- gseKEGG(geneList=geneList, organism='hsa', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE, scoreType="pos")
                if (!is.null(gse_kegg_results) && nrow(as.data.frame(gse_kegg_results)) > 0) {
                    p4 <- dotplot(gse_kegg_results, showCategory=15)
                    ggsave(file.path(output_dir, "plot_kegg_dotplot.png"), plot = p4, width=10, height=8)
                    p5 <- ridgeplot(gse_kegg_results, showCategory=15)
                    ggsave(file.path(output_dir, "plot_kegg_ridgeplot.png"), plot = p5, width=10, height=8)
                    p6 <- gseaplot2(gse_kegg_results, geneSetID = 1:min(3, nrow(as.data.frame(gse_kegg_results))))
                    ggsave(file.path(output_dir, "plot_kegg_gseaplot.png"), plot = p6, width=10, height=8)
                }
            }, error = function(e) {
                print(paste("An error occurred during GSEA analysis:", e$message))
            })
        }
        """

        # Step 3: R 코드 실행
        st.info("Python 내에서 R 코드를 직접 실행하여 GSEA 분석을 시작합니다...")
        robjects.r(r_code)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            robjects.r['run_gsea_in_r'](r_df, temp_dir)
            st.success("GSEA 분석이 성공적으로 완료되었습니다!")

            plots = {}
            plot_files = ["plot_go_dotplot.png", "plot_go_ridgeplot.png", "plot_go_gseaplot.png", 
                          "plot_kegg_dotplot.png", "plot_kegg_ridgeplot.png", "plot_kegg_gseaplot.png"]
            for plot_name in plot_files:
                plot_path = os.path.join(temp_dir, plot_name)
                if os.path.exists(plot_path):
                    plots[plot_name] = plot_path
            return plots
            
    except Exception as e:
        st.error(f"rpy2를 이용한 R 코드 실행 중 오류가 발생했습니다: {e}")
        return None

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
                
                plots = process_and_run_gsea_rpy2(prescription_name_input, selected_herbs_info, herb_weights)
                
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

