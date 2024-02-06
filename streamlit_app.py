import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import math


st.set_page_config(
    page_title="Gene Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://www.extremelycoolapp.com/help',
        'Report a bug': "https://www.extremelycoolapp.com/bug",
        'About': "# This is a header. This is an *extremely* cool app!"
    }
)

st.title("Spatial Transcriptomics Analysis of P3")
st.header("Gene Dashboard")
st.subheader("Perony Nogueira, PhD Candidate")
leaflets_data = pd.read_csv("leaflets_data.csv")

# To find the gene of interest (should be precise)
def filter_genes(geneName):
  gene_entries = leaflets_data[leaflets_data['Gene'] == geneName]
  return gene_entries


gene_name_list = list(set(leaflets_data["Gene"]))
gene_name_list.insert(0,"")

tab1, tab2 = st.tabs(["Genes Comparison", "Top Correlation"])

with tab1:
    with st.form("Genes Selection",clear_on_submit=True):
        gName1 = st.selectbox("Select 1st Gene", options=gene_name_list)
        gName2 = st.selectbox("Select 2nd Gene", options=gene_name_list)
        submitted = st.form_submit_button("Submit")
        if submitted:

            gene1 = filter_genes(gName1)
            gene2 = filter_genes(gName2)
            # Extract the relevant columns from the data frame
            gene1_data = gene1[['BarcodeID', 'Read Count']]
            gene2_data = gene2[['BarcodeID', 'Read Count']]

            # Merge the two data frames on the 'BarcodeID' column
            merged_data = pd.merge(gene1_data, gene2_data, on='BarcodeID')

            # Calculate the Spearman correlation coefficient between the read counts of the two genes
            correlation, pvalue = spearmanr(merged_data['Read Count_x'], merged_data['Read Count_y'])

            # Print the correlation coefficient
            print('Spearman correlation coefficient:', correlation)
            print("p-value:", pvalue)



            # Assuming 'merged_data' is your DataFrame containing 'Read Count_x' and 'Read Count_y'
            fig = px.scatter(merged_data, x='Read Count_x', y='Read Count_y', labels={'Read Count_x': gName1, 'Read Count_y': gName2})

            textbox_str = f"Correlation: {correlation:.4f}<br>P-value: {pvalue}"

            fig.add_annotation(
                go.layout.Annotation(
                    text=textbox_str,
                    align='left',
                    showarrow=False,
                    xref='paper',
                    yref='paper',
                    x=0.05,
                    y=0.95,
                    bordercolor='black',
                    borderwidth=1,
                    bgcolor='#dddddd',
                    opacity=0.9,
                    font=dict(
                        color='#141414'
                    )
                )
            )

            # Customize layout as needed
            fig.update_layout(
                title=f'Correlation between {gName1} and {gName2}',
                xaxis_title=gName1,
                yaxis_title=gName2
            )

            # Show the plot
            st.plotly_chart(fig)

with tab2:
    TOPcorrelations = dict()
    with st.form("Top Correlations",clear_on_submit=False):
        gName1 = st.selectbox("Select the gene", options=gene_name_list)
        correlation_threshold = st.slider("Correlation Threshold",0.2,0.8,0.6)
        pvalue_threshold = st.slider("P-value Threshold",0.001,0.05,0.01)
        submit2 = st.form_submit_button("Analyze")
        if submit2:
            gene1 = filter_genes(gName1)
            gene1_data = gene1[['BarcodeID', 'Read Count']]
            with st.spinner("🔬 Analyzing Gene Correlations, hang tight... ⏳"):
                for i in gene_name_list[1:]:
                    if i != gName1:
                        # st.write("Gene", i)
                        gene2 = filter_genes(i)
                        gene2_data = gene2[['BarcodeID', 'Read Count']]
                        merged_data = pd.merge(gene1_data, gene2_data, on='BarcodeID')
                        correlation, pvalue = spearmanr(merged_data['Read Count_x'], merged_data['Read Count_y'])
                        if (not math.isnan(correlation)) and (not math.isnan(pvalue)):
                            if (correlation >= correlation_threshold) and pvalue <= pvalue_threshold:
                                key_name = f"{gName1} and {i}"
                                TOPcorrelations[key_name] = (correlation,pvalue)
                # st.write(TOPcorrelations)

                # Convert dictionary to DataFrame
                df_top = pd.DataFrame(list(TOPcorrelations.items()), columns=['Gene', 'Values'])

                # Split the 'Values' tuple into two separate columns
                df_top[['Correlation', 'PValue']] = pd.DataFrame(df_top['Values'].tolist(), index=df_top.index)

                # Drop the 'Values' column as it's no longer needed
                df_top.drop(columns=['Values'], inplace=True)
                st.dataframe(df_top)
                st.balloons()
                st.success("✅ Processing Complete! All correlations have been successfully identified. Thank you for your patience. 🎉")
                # df_top.to_csv(f"TOP_Correlations_with_{gName1}.csv")