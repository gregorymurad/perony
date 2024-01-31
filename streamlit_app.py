import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

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
                bgcolor='lightgray',
                opacity=0.9
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