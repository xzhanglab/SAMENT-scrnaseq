<p align="center">
  <img src="GitHub.png" alt="Project banner" style="width:100%; max-width:100%; height:auto;" />
</p>

# Unbiased niche labeling maps immune-excluded niche in bone metastasis

This repository is the manuscript companion resource for the study above. It collects the analysis notebooks, R Markdown workflows, interactive pathway viewers, and Zenodo-hosted processed datasets used to reproduce the single-cell, pathway, and supporting analyses from the paper.

## Companion spatial analysis repository

The image-based spatial analysis used for this study has now been organized as a separate companion application:

- **Companion app:** [TME Spatial](https://github.com/fengshuoliu/TME_spatial)
- **Companion website:** [https://fengshuoliu.github.io/TME_spatial/](https://fengshuoliu.github.io/TME_spatial/)

### How the two repositories connect

| Resource | Primary purpose | Use this when you want to... |
| --- | --- | --- |
| **This repository** | Reproduce the manuscript-level single-cell, pathway, trajectory, and communication analyses | Follow the paper workflow, inspect the scripts behind the figures, and download the manuscript companion datasets |
| **TME Spatial** | Run the immunofluorescence image spatial-analysis workflow as a shared macOS / Windows app | Re-run or adapt the image-based spatial analysis, niche labeling, ROI analysis, cell distribution analysis, and distance analysis used in this study |

If you are mainly interested in the **image-based spatial analysis** from this paper, start with **TME Spatial**.  
If you want the broader **paper reproduction resource**, including scRNA-seq and pathway workflows, start with **this repository**.

## Project links

| Link | Description |
| --- | --- |
| [Project website](https://xzhanglab.github.io/SAMENT-scrnaseq/) | Static website for this manuscript resource |
| [Zenodo dataset collection](https://doi.org/10.5281/zenodo.14796581) | Downloadable processed data and demonstration archives |
| [TME Spatial repository](https://github.com/fengshuoliu/TME_spatial) | Companion repository for the image-based spatial-analysis app |
| [TME Spatial website](https://fengshuoliu.github.io/TME_spatial/) | Installation, workflow, and documentation for the companion app |

## Citation

If you use this resource or the companion spatial-analysis workflow, please cite:

Xu Z, Liu F, Ding Y, Pan T, Wu Y-H, Han Y, Liu J, Bado IL, Zhang W, Wu L, Gao Y, Hao X, Yu L, Li Y, Edwards DG, Chan HL, Aguirre S, Dieffenbach MW, Chen E, Wang S, Shen Y, Hoffman D, Becerra Dominguez L, Rivas CH, Chen X, Wang H, Kang Y, Gugala Z, Satcher RL, Zhang XH-F. *Unbiased niche labeling maps immune-excluded niche in bone metastasis.* Cell. 2026. Published online April 2026. doi:10.1016/j.cell.2026.04.009

## Interactive resources

### Pathway exploration apps

Interactive volcano plots allow exploration of differentially regulated pathways between biotin-positive and biotin-negative populations.

- **SAMENT macrophage explorer:** https://samentexplore-73cfuwnsd8tzxatzwvcwvv.streamlit.app
- **SAMENT neutrophil explorer:** https://neutrophilbiotinpositive-vs-negativegsvapy-dd42nk8mhnrf4vjhvuo.streamlit.app/

<p align="center">
  <img src="SAMENT_explore/demo.gif" width="900" alt="SAMENT explorer demo" />
</p>

### Spatial image analysis app

The companion **TME Spatial** app was developed from the spatial-analysis needs of this study and now provides a dedicated workflow for:

- multi-channel image input and configuration
- overlay and split-channel visualization
- nuclei segmentation
- cell type assignment
- neighborhood clustering
- ROI / region analysis
- cell distribution analysis
- nearest-neighbor and cell-to-boundary distance analysis

Use TME Spatial together with the Zenodo spatial-analysis demonstration archives listed below.

## Overview

This repository guides you through the following:

1. Integrating datasets, applying batch correction, and reproducing the single-cell analyses from the manuscript.
2. Reproducing the integrated bulk / microarray and pathway analyses.
3. Accessing the companion image-based spatial-analysis resources and connecting them to the dedicated TME Spatial app.

## Analysis pipeline files

| File name | Description |
| --- | --- |
| `01_preprocessing_demultiplexing_template.Rmd` | Process Cell Ranger outputs and generate individual Seurat objects, including demultiplexing |
| `02_integration.Rmd` | Integrate the single-cell objects |
| `03_scanpy_plot.ipynb` | Generate manuscript plots in Scanpy |
| `04_processing_for_DESeq2.Rmd` | Prepare DESeq2 objects for differential-expression and pathway analysis |
| `05_DESeq2_DEG.Rmd` | Differential gene-expression analysis |
| `06.1_pre_requested_Matrix.utils.R` | Utility script for required packages and helper functions |
| `06.2_GSEA.Rmd` | Gene set enrichment analysis |
| `07.GSVA.Rmd` | GSVA pathway activity analysis |
| `08.velocyto.ipynb` | RNA velocity / trajectory analysis |
| `09.expression_distence.ipynb` | Transcriptomic distance analysis |
| `10.CellChat_comparsion.Rmd` | Cell-cell communication analysis |

## Data files from Zenodo

All intermediate data produced by the workflows in this project are available on [Zenodo](https://doi.org/10.5281/zenodo.14796581).

The files most directly connected to **TME Spatial** are the immunofluorescence image spatial-analysis archives.

| Directory or file | Description | Related figure(s) | Connection to TME Spatial |
| --- | --- | --- | --- |
| `immunofluourscance-image-spatial-analysis_ERE_or_ERa-positive-macrophages_distribution.zip` | Spatial analysis measuring ERa+ / ERE+ macrophage distance to nearest tumor cells | Fig S4G-L | Can be re-run or adapted in the companion spatial-analysis workflow |
| `immunofluourscance-image-spatial-analysis_diffusion-lesions_niche-labeling-model-comparisions.zip` | Spatial analysis comparing niche-labeling methods in bone metastasis diffusion lesions | Fig 2E, Fig S2B, Fig S2C | Closely related to the image-based niche-labeling workflow in TME Spatial |
| `immunofluourscance_image_spatial_analysis_demonstration_scripts_and_data.zip` | Demonstration scripts and data for customized spatial analysis | / | Best entry point for connecting this paper resource to TME Spatial |
| `immunofluourscance-image-spatial-analysis_isolate-lesions_niche-labeling-model-comparisions.zip` | Spatial analysis comparing niche-labeling methods in isolated lesions | Fig 2D, Fig S2A | Related to the image-based niche-labeling workflow in TME Spatial |
| `immunofluourscance-image-spatial-analysis_T_cell_distribution_in_ctrl-Esr1KO.zip` | Spatial analysis measuring T cell distribution in the tumor microenvironment | Fig 7 | Related to cell distribution analysis in TME Spatial |
| `SAMENT_single_cell_per-sample-Seurat_objects.zip` | Individual SAMENT Seurat and Scanpy objects from demultiplexed scRNA-seq, without cell-type annotation | / | This repository |
| `SAMENT_single_cell_integrated_objects.zip` | Integrated SAMENT Seurat and Scanpy objects, batch-corrected and cell-type annotated | All SAMENT scRNA-seq related panels | This repository |
| `cellranger_demultiplexed_scRNA-seq_matrix.zip` | Cell Ranger scRNA-seq outputs | / | This repository |
| `Macrophage_LysM_Esr1_KO_single_cell_integrated_objects.zip` | Integrated LysM-Esr1 Seurat and Scanpy objects, batch-corrected and cell-type annotated | All LysM-Esr1 scRNA-seq related panels | This repository |
| `Macrophage_LysM_Esr1_KO_single_cell_per-sample-Seurat_objects.zip` | Individual LysM-Esr1 Seurat and Scanpy objects from demultiplexed scRNA-seq, without cell-type annotation | / | This repository |
| `scRNA-seq_analysis_scripts.zip` | All scRNA-seq analysis scripts | / | This repository |

## Suggested wording for linking the two repositories

If you want the two repositories to feel explicitly connected, this short paragraph works well in both places:

> This manuscript resource provides the paper-level analysis workflows and processed data, while the companion repository **TME Spatial** packages the image-based spatial-analysis workflow used in this study into a reusable application for macOS and Windows.

## Prepared by

Fengshuo Liu
