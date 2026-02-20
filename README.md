<!-- Top banner image (full width) -->
<p align="center">
  <img src="GitHub.png" alt="Project banner" style="width:100%; max-width:100%; height:auto;" />
</p>


# <span style="font-size: 16px;">Unbiased metastatic niche-labeling identifies estrogen receptor-positive macrophages as a barrier of T cell infiltration during bone colonization</span>

<span style="font-size: 12px;">
This repository provides instructions and code to reproduce the major results, numerics, and figures from the <a href="https://doi.org/10.1101/2024.05.07.593016"><b>manuscript</b></a>:
</span>

### <span style="font-size: 14px;">Citation</span>
<span style="font-size: 12px;">
Xu, Z., Liu, F., Ding, Y., Hao, X., Pan, T., Miles, G., Wu, Y.-H., Liu, J., Bado, I. L., Zhang, W., Wu, L., Gao, Y., Yu, L., Edwards, D. G., Chan, H. L., Aguirre, S., Dieffenbach, M. W., Chen, E., Shen, Y., Hoffman, D., Dominguez, L. B., Rivas, C. H., Chen, X., Wang, H., Gugala, Z., Satcher, R. L., & Zhang, X. H.-F. (2024). Unbiased metastatic niche-labeling identifies estrogen receptor-positive macrophages as a barrier of T cell infiltration during bone colonization. <i>Cell</i> (accepted). https://doi.org/10.1101/2024.05.07.593016
</span>


<br>

### <span style="font-size: 14px;">Interactive Volcanoplot: exploying differentially regulated pathways between biotin positive and negative populations</span>
#### It may take a while to load the app, please be patient
<span style="font-size: 12px;">
SAMENT Macrophage: https://samentexplore-73cfuwnsd8tzxatzwvcwvv.streamlit.app
<br>

<span style="font-size: 12px;">
SAMENT Neutrophil: https://neutrophilbiotinpositive-vs-negativegsvapy-dd42nk8mhnrf4vjhvuo.streamlit.app/
</a>
<br><br>
<img src="SAMENT_explore/demo.gif" width="800" alt="SAMENT Neutrophil Demo"/>
</span>
  
### <span style="font-size: 14px;">Data</span>
<span style="font-size: 12px;">
All intermediate data produced by running this code, as described below, are available for download on <a href="https://doi.org/10.5281/zenodo.14796581"><b>Zenodo</b></a>.
</span>

### <span style="font-size: 14px;">Overview</span>
<span style="font-size: 12px;">
These instructions will guide you through the following:
<br>
1. Integrating datasets, applying batch correction, and reproducing analysis from the manuscript.  
<br>
2. Reproducing the results from analyzing integrated bulk/microarray datasets.  
<br>
</span>


| **File Name**                                    | **Description**                                                                 | 
|--------------------------------------------------|---------------------------------------------------------------------------------|
| <span style="font-size: 12px;">01_preprocessing_demultiplexing_template.Rmd</span>      | <span style="font-size: 12px;">Process Cell Ranger outputs and generate individual Seurat objects, including demultiplexing</span>              |
| <span style="font-size: 12px;">02_integration.Rmd</span>          | <span style="font-size: 12px;">Integration of the single cell objects</span>                  |
| <span style="font-size: 12px;">03_scanpy_plot.ipynb</span>                                    | <span style="font-size: 12px;">Generate the plots in Scanpy</span>                                              | 
| <span style="font-size: 12px;">04_processing_for_DESeq2.Rmd</span>     | <span style="font-size: 12px;">Prepare DESeq2 object for DEG and pathway analysis</span>                              |
| <span style="font-size: 12px;">05_DESeq2_DEG.Rmd</span>                       | <span style="font-size: 12px;">DEG analysis</span>                                                                   |         |
| <span style="font-size: 12px;">06.1_pre_requested_Matrix.utils.R</span>                                      | <span style="font-size: 12px;">Prerequested package</span>                                                   |                 |
| <span style="font-size: 12px;">06.2_GSEA.Rmd</span>                                      | <span style="font-size: 12px;">GSEA analysis</span>                                                                  |
| <span style="font-size: 12px;">07.GSVA.Rmd</span>                       | <span style="font-size: 12px;">GSVA analysis</span>                                                           | 
| <span style="font-size: 12px;">08.velocyto.ipynb</span>                                  | <span style="font-size: 12px;">trajectory analysis</span>                                               |
| <span style="font-size: 12px;">09.expression_distence.ipynb</span> | <span style="font-size: 12px;">Transcriptomic differences</span>                                   |
| <span style="font-size: 12px;">10.CellChat_comparsion.Rmd</span> | <span style="font-size: 12px;">Cell-cell conmunication analysis</span>                                   |

---
### <span style="font-size: 14px;">Data files (from Zenodo)</span>

<style>
/* Responsive table for GitHub README */
.zenodo-table {
  width: 100%;
  table-layout: fixed;          /* key: respect column widths */
  border-collapse: collapse;
}

.zenodo-table th,
.zenodo-table td {
  padding: 6px 10px;
  vertical-align: top;
  font-size: 12px;
}

/* column widths (adjust as you like) */
.zenodo-table col:nth-child(1) { width: 30%; }
.zenodo-table col:nth-child(2) { width: 55%; }
.zenodo-table col:nth-child(3) { width: 15%; }

/* wrap long filenames nicely */
.zenodo-table td:first-child,
.zenodo-table th:first-child {
  overflow-wrap: anywhere;      /* break long tokens */
  word-break: break-word;
}

/* optional: keep Related Figure from wrapping too much */
.zenodo-table td:last-child,
.zenodo-table th:last-child {
  white-space: nowrap;
}
</style>

<table class="zenodo-table">
  <colgroup>
    <col><col><col>
  </colgroup>
  <thead>
    <tr>
      <th><b>Directory/File</b></th>
      <th><b>Description</b></th>
      <th><b>Related Figure</b></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>immunofluourscance-image-spatial-analysis_ERE_or_ERa-positive-macrophages_distribution.zip</td>
      <td>Spatial analysis measures Erα+/ERE+ macrophage distance to its nearest tumor cell.</td>
      <td>Fig S4G–L</td>
    </tr>
    <tr>
      <td>immunofluourscance-image-spatial-analysis_diffusion-lesions_niche-labeling-model-comparisions.zip</td>
      <td>Spatial analysis compares niche-labeling methods in bone metastasis (diffusion lesions).</td>
      <td>Fig 2E, Fig S2B, Fig S2C</td>
    </tr>
    <tr>
      <td>SAMENT_single_cell_per-sample-Seurat_objects.zip</td>
      <td>Individual SAMENT Seurat and Scanpy objects from demultiplexed scRNA-seq (no cell-type annotation).</td>
      <td>/</td>
    </tr>
    <tr>
      <td>immunofluourscance_image_spatial_analysis_demonstration_scripts_and_data.zip</td>
      <td>Demonstration scripts and data for customized spatial analysis.</td>
      <td>/</td>
    </tr>
    <tr>
      <td>SAMENT_single_cell_integrated_objects.zip</td>
      <td>Integrated SAMENT Seurat and Scanpy objects (batch-corrected, cell-type annotated).</td>
      <td>All SAMENT scRNA-seq related panels</td>
    </tr>
    <tr>
      <td>immunofluourscance-image-spatial-analysis_isolate-lesions_niche-labeling-model-comparisions.zip</td>
      <td>Spatial analysis compares niche-labeling methods in bone metastasis (isolated lesions).</td>
      <td>Fig 2D, Fig S2A</td>
    </tr>
    <tr>
      <td>cellranger_demultiplexed_scRNA-seq_matrix.zip</td>
      <td>scRNA-seq Cell Ranger outputs.</td>
      <td>/</td>
    </tr>
    <tr>
      <td>Macrophage_LysM_Esr1_KO_single_cell_integrated_objects.zip</td>
      <td>Integrated LysM-Esr1 Seurat and Scanpy objects (batch-corrected, cell-type annotated).</td>
      <td>All LysM-Esr1 scRNA-seq related panels</td>
    </tr>
    <tr>
      <td>Macrophage_LysM_Esr1_KO_single_cell_per-sample-Seurat_objects.zip</td>
      <td>Individual LysM-Esr1 Seurat and Scanpy objects from demultiplexed scRNA-seq (no cell-type annotation).</td>
      <td>/</td>
    </tr>
    <tr>
      <td>immunofluourscance-image-spatial-analysis_T_cell_distribution_in_ctrl-Esr1KO.zip</td>
      <td>Spatial analysis measures T cell distribution in the TME.</td>
      <td>Fig 7</td>
    </tr>
    <tr>
      <td>scRNA-seq_analysis_scripts.zip</td>
      <td>All scRNA-seq analysis scripts.</td>
      <td>/</td>
    </tr>
  </tbody>
</table>


---


### <span style="font-size: 14px;">Prepared by Fengshuo Liu</span>
