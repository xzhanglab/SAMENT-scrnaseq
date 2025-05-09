# <span style="font-size: 16px;">Unbiased metastatic niche-labeling identifies estrogen receptor-positive macrophages as a barrier of T cell infiltration during bone colonization</span>

<span style="font-size: 12px;">
This repository provides instructions and code to reproduce the major results, numerics, and figures from the <a href="https://doi.org/10.1101/2024.05.07.593016"><b>manuscript</b></a>:
</span>

### <span style="font-size: 14px;">Citation</span>
<span style="font-size: 12px;">
Liu, F., Xu, Z., Ding, Y., Hao, X., Pan, T., Miles, G., Wu, Y.-H., Liu, J., Bado, I. L., Zhang, W., Wu, L., Gao, Y., Yu, L., Edwards, D. G., Chan, H. L., Aguirre, S., Dieffenbach, M. W., Chen, E., Shen, Y., Hoffman, D., Dominguez, L. B., Rivas, C. H., Chen, X., Wang, H., Gugala, Z., Satcher, R. L., & Zhang, X. H.-F. (2024).Unbiased metastatic niche-labeling identifies estrogen receptor-positive macrophages as a barrier of T cell infiltration during bone colonization. bioRxiv. https://doi.org/10.1101/2024.05.07.593016
</span>

<br>

### <span style="font-size: 14px;">Inter-active Volcanoplot: exploying differentially regulated pathways between biotin positive and negative populations</span>
#### it may take a while to load the app, please be patient
<span style="font-size: 12px;">
SAMENT Macrophage: https://samentexplore-jxgg9fba6fuyn38wpxvwc7.streamlit.app
<br>

<span style="font-size: 12px;">
SAMENT Neutrophil: https://neutrophilbiotinpositive-vs-negativegsvapy-dd42nk8mhnrf4vjhvuo.streamlit.app/
<br>
  
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

## <span style="font-size: 14px;">Data file (from Zenodo)</span>

| Directory/File                      | Description                                                                 |
|-------------------------------------|-----------------------------------------------------------------------------|
| <span style="font-size: 12px;">analysis_scripts.zip</span>     | <span style="font-size: 12px;">Zenodo deposit of the analysis script, same files in GitHub</span>             |
| <span style="font-size: 12px;">SAMENT_single_cell_integrated_objects.zip</span>         | <span style="font-size: 12px;">Integrated single cell objects from SAMENT scRNA-seq, including both Seurat and Scanpy formats</span>                    |
| <span style="font-size: 12px;">SAMENT_single_cell_per-sample-Seurat_objects.zip</span>                              | <span style="font-size: 12px;">Individual per-sample Seurat objects from SAMENT scRNA-seq</span> |
| <span style="font-size: 12px;">Macrophage_LysM_Esr1_KO_single_cell_integrated_objects.zip</span>           | <span style="font-size: 12px;">Integrated single cell objects from scRNA-seq for in vivo macrophage-specific Esr1 knock down, including both Seurat and Scanpy formats</span> |
| <span style="font-size: 12px;">Macrophage_LysM_Esr1_KO_single_cell_per-sample-Seurat_objects.zip</span>                         | <span style="font-size: 12px;">Individual per-sample Seurat objects from scRNA-seq, for macrophage-specific Esr1 knock out mice</span>              |
| <span style="font-size: 12px;">Macrophage_LysM_Esr1_KO_CellChat_objects.zip</span>     | <span style="font-size: 12px;">Cell-Cell conmunication objects</span> |

---

### <span style="font-size: 14px;">Prepared by Fengshuo Liu</span>
