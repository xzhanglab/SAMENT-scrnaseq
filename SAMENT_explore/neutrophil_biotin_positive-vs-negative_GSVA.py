#!/usr/bin/env python
# coding: utf-8

# In[5]:
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
from plotly.io import to_image
import textwrap

# Check if kaleido is installed for image export
try:
    import kaleido
    kaleido_available = True
except ImportError:
    kaleido_available = False

# Load the data (Cached for performance)
@st.cache_data
def load_data():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'macrophages_biotin_positive-vs-negative_GSVA.csv')
    
    if not os.path.exists(file_path):
        st.error(f"File {file_path} not found. Please check the file path.")
        return None
    
    df = pd.read_csv(file_path)
    df = df.set_index(df.columns[0])
    df.index = df.index.str.strip()
    df['-log10(adj.P.Val)'] = -np.log10(df['P.Value'])
    return df

df = load_data()

# MSigDB table for Prefix selection
msigdb_table = pd.DataFrame({
    'Category': ['H: Hallmark Gene Sets', 'C1: Positional Gene Sets', 'C2:CP (Canonical Pathways)', 
                 'C2:CGP (Chemical and Genetic Perturbations)', 'C3:MIR (MicroRNA targets)', 'C3:TFT (Transcription Factor targets)',
                 'C4:CGN (Cancer Gene Neighborhoods)', 'C4:CM (Cancer Modules)', 'C5: GOBP (Biological Processes)',
                 'C5: GOCC (Cellular Components)', 'C5: GOMF (Molecular Functions)', 'C6: Oncogenic Signatures',
                 'C7: Immunologic Signatures', 'C8: Cell Type Signatures'],
    'Prefix': ['HALLMARK', 'CHR', 'KEGG, REACTOME, BIOCARTA', 'CGP', 'MIR', 'TFT',
               'CGN', 'CM', 'GOBP', 'GOCC', 'GOMF', 'C6, ONCOGENIC', 'C7', 'C8'],
    'Example': ['HALLMARK_APOPTOSIS', 'CHR1Q22', 'KEGG_APOPTOSIS, REACTOME_CELL_CYCLE', 'CGP_CANCER_DRUGS', 'MIR_21', 'TFT_STAT3',
                'CGN_P53', 'CM_APOPTOSIS', 'GOBP_APOPTOSIS', 'GOCC_NUCLEUS', 'GOMF_RECEPTOR_ACTIVITY', 'C6_MYC_TARGETS',
                'C7_T_CELL_RECEPTOR_PATHWAY', 'C8_T_CELL_SIGNATURE']
})

# Function to categorize pathways
def get_category(row, keywords=[], exclude_keywords=[], logic='AND'):
    pathway_name = row.name.replace('_', ' ').upper()
    keywords = [kw.upper().strip() for kw in keywords if kw.strip() != '']
    exclude_keywords = [kw.upper().strip() for kw in exclude_keywords if kw.strip() != '']
    
    # Exclude pathways that match any excluded keywords
    if any(keyword in pathway_name for keyword in exclude_keywords):
        return 'non-significant'
    
    if logic == 'AND':
        if all(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    elif logic == 'OR':
        if any(keyword in pathway_name for keyword in keywords):
            return 'keyword_match'
    
    if row['GSVA_score'] > 0.2 and row['-log10(adj.P.Val)'] > 1:
        return 'upregulated'
    elif row['GSVA_score'] < -0.2 and row['-log10(adj.P.Val)'] > 1:
        return 'downregulated'
    else:
        return 'non-significant'

# Function to wrap long text in legend
def wrap_text(text, width=30):
    return '<br>'.join(textwrap.wrap(text, width=width))

# Function to update the plot
def update_plot(keywords=[], exclude_keywords=[], logic='AND', width='100%', height=800, interactive=True):
    df['category'] = df.apply(get_category, axis=1, keywords=keywords, exclude_keywords=exclude_keywords, logic=logic)
    palette = {'keyword_match': '#32CD32', 'upregulated': '#FF6347', 'downregulated': '#1E90FF', 'non-significant': '#A9A9A9'}
    fig = go.Figure()
    
    # Plot non-significant pathways
    non_significant_df = df[df['category'] == 'non-significant']
    fig.add_trace(go.Scatter(x=non_significant_df['GSVA_score'], y=non_significant_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['non-significant'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[name for name in non_significant_df.index], hoverinfo='text', name='Non-Significant'))
    
    # Plot up-regulated pathways
    upregulated_df = df[df['category'] == 'upregulated']
    fig.add_trace(go.Scatter(x=upregulated_df['GSVA_score'], y=upregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['upregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["upregulated"]};">{name}</span>' for name in upregulated_df.index], hoverinfo='text', name='Up-regulated'))

    # Plot down-regulated pathways
    downregulated_df = df[df['category'] == 'downregulated']
    fig.add_trace(go.Scatter(x=downregulated_df['GSVA_score'], y=downregulated_df['-log10(adj.P.Val)'], mode='markers',
                             marker=dict(size=8, color=palette['downregulated'], opacity=0.8, line=dict(width=0.5, color='black')),
                             text=[f'<span style="color:{palette["downregulated"]};">{name}</span>' for name in downregulated_df.index], hoverinfo='text', name='Down-regulated'))
    
    # Sort keyword-matched pathways by P.Value
    keyword_df = df[df['category'] == 'keyword_match'].sort_values('P.Value')
    
    if interactive:
        # Plot keyword-matched pathways interactively
        fig.add_trace(go.Scatter(x=keyword_df['GSVA_score'], y=keyword_df['-log10(adj.P.Val)'], mode='markers',
                                 marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                 text=[f'<span style="color:{palette["keyword_match"]};">{name}</span>' for name in keyword_df.index],
                                 hoverinfo='text', name=wrap_text(', '.join(keywords), width=20)))
    else:
        # Plot numbered keyword-matched pathways (Start numbering from 1)
        for i, (index, row) in enumerate(keyword_df.iterrows(), start=1):
            showlegend = i == 1  # Only show legend for the first trace
            fig.add_trace(go.Scatter(x=[row['GSVA_score']], y=[row['-log10(adj.P.Val)']], mode='text+markers',
                                     marker=dict(size=15, color=palette['keyword_match'], opacity=0.8, line=dict(width=0.5, color='black')),
                                     text=f"<b style='color:black;'>{i}</b>",  # Bold and black color for numbers, starting from 1
                                     hoverinfo='text', name=wrap_text(f"{', '.join(keywords)}") if showlegend else None, showlegend=showlegend))

    # Add vertical dashed lines at x = -0.2 and x = 0.2
    fig.update_layout(
        shapes=[
            dict(
                type="line",
                x0=-0.2,
                y0=0,
                x1=-0.2,
                y1=df['-log10(adj.P.Val)'].max(),
                line=dict(
                    color="grey",
                    width=2,
                    dash="dash",
                ),
                layer='above'  # Ensure the line is drawn above the data points
            ),
            dict(
                type="line",
                x0=0.2,
                y0=0,
                x1=0.2,
                y1=df['-log10(adj.P.Val)'].max(),
                line=dict(
                    color="grey",
                    width=2,
                    dash="dash",
                ),
                layer='above'
            )
        ]
    )

    # Set layout with wrapped legend text
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(255,255,255,1)',
        title='Macrophage: Biotin Positive vs Negative',
        xaxis_title='GSVA Score',
        yaxis_title='-log10(adj.P.Val)',
        title_font_size=18,
        width=width,
        height=height,
        legend_title_text='Pathway Categories',
        autosize=True,
        legend=dict(
            x=1, y=1, 
            bgcolor='rgba(255,255,255,0)',  # Transparent background for legend
            bordercolor='rgba(255,255,255,0)',  # Transparent border
        )
    )
    
    return fig, keyword_df

if df is not None:
    # Sidebar input for keywords and logic
    st.sidebar.header('Search Parameters')
    num_keywords = st.sidebar.number_input('Number of Keywords', min_value=1, max_value=10, value=2)
    keywords = [st.sidebar.text_input(f'Keyword {i+1}') for i in range(num_keywords)]
    logic = st.sidebar.selectbox('Logic', ['AND', 'OR'])
    exclude_keywords = st.sidebar.text_area('Exclude Keywords (comma-separated)').split(',')
    keywords = [kw for kw in keywords if kw.strip() != '']
    
    fig_width = st.sidebar.slider('Figure Width', min_value=400, max_value=1600, value=1000, step=50)
    fig_height = st.sidebar.slider('Figure Height', min_value=400, max_value=1200, value=800, step=50)
    
    # Prefix selection for filtering pathways
    st.sidebar.header('Filter Pathways by Prefix')
    selected_prefixes = st.sidebar.multiselect('Select Prefixes', msigdb_table['Prefix'].tolist(), default=msigdb_table['Prefix'].tolist())
    
    # Filter the data based on selected prefixes
    if selected_prefixes:
        pattern = '|'.join(selected_prefixes)  # Create a pattern to match selected prefixes
        df = df[df.index.str.contains(pattern)]

    interactive_keywords = st.sidebar.radio('Keyword-Matched Pathways Interactive?', ('Yes', 'No'))
    
    # Update plot based on user input
    fig, keyword_df = update_plot(keywords, exclude_keywords, logic, width=fig_width, height=fig_height, interactive=(interactive_keywords == 'Yes'))
    
    # Display plot
    st.plotly_chart(fig, use_container_width=True)

    # Display the table for keyword-matched pathways if interactive is set to No
    if interactive_keywords == 'No' and not keyword_df.empty:
        st.write("### Keyword-Matched Pathways")
        keyword_df_display = keyword_df[['P.Value']].reset_index().rename(columns={'index': 'Pathway'})
        keyword_df_display.index += 1  # Ensure the table starts numbering from 1
        st.dataframe(keyword_df_display)

    # Display the MSigDB Table
    st.write("### MSigDB Categories And Prefixes")
    st.table(msigdb_table)
    st.write("You can specify the category (\"Prefix\") desired for your research in search box.")

    
    # Download plot as PNG or PDF
    st.sidebar.header('Download Plot')
    download_format = st.sidebar.radio('Download Format', ('PNG', 'PDF'))
    
    if st.sidebar.button('Download'):
        if kaleido_available:
            if download_format == 'PNG':
                file_bytes = to_image(fig, format='png', engine="kaleido", scale=3)  # 300 DPI
                st.sidebar.download_button(label='Download as PNG', data=file_bytes, file_name='plot.png', mime='image/png')
            elif download_format == 'PDF':
                file_bytes = to_image(fig, format='pdf', engine="kaleido", scale=3)  # 300 DPI
                st.sidebar.download_button(label='Download as PDF', data=file_bytes, file_name='plot.pdf', mime='application/pdf')
        else:
            st.sidebar.error("Image export requires the 'kaleido' package. Please install it by adding 'kaleido' to your requirements.txt file.")
