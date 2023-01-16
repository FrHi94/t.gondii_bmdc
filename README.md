# scDual-Seq of Toxoplasma gondii-infected mouse bone marrow-derived dendritic cells reveals host cell heterogeneity and differential infection dynamics

## Summary 

Host-parasite interactions include complex interplays between an invading and a defending organism, each continuously adapting to gain the upper hand. The protozoan parasite, Toxoplasma gondii, can invade every nucleated cell type in a vertebrate host, including immune cells while the host mounts a protective response. 
Here, we utilize Dual-scSeq to parse out heterogeneous transcription of bone marrow-derived dendritic cells (BMDCs) infected with T. gondii type I, RH (LDM) or type II, ME49 (PTG) parasites, over multiple time points post infection.
We find that the two parasite lineages distinctly manipulate two subpopulations of infected BMDCs. Co-expression networks establish host and parasite genes, with implications for modulation of host immunity and host-pathogen interactions. Integration of published data validates immune pathways and suggests novel candidate genes involved in host-pathogen interactions. This study aims to provide a comprehensive resource for future characterization of host-pathogen interplay among other protozoan parasites within their host niches, as well as that of bacterial and viral pathogens. 


## Structure 

**`data`** - contains count matrices and metadata 

   - **`FACS_sorting`** contains files describing gating strategies for sorting of infected BMDCs in '.docx' format

   - **`cnt`** contains

   - raw (non-normalized) count-matrix after mapping ('cnt_mus_toxo_raw.csv') after filtering cells

   - count-matrix from microarray experiments ('GSE134869_SP7_umi_counts.tsv.gz') by Pandey et al. `(1)`

   - normalized read counts for mouse data ('mouse-DC-counts-filter.csv')

   - normalized read counts for *Toxoplasma gondii* data ('toxo-DC-counts.csv')


   - **`cnt/feat_count`** contains raw (non-normalized) count-matrices for individual plates used for smartSeq2 protocol


**`res`** - contains intermediate results from the R markdown scripts needed for some of the analyses 

   - **`GEO_comparison`** contains `.tsv`  files of differentially expressed genes by all single conditions tested in Pandey et al. study `(1)`

   - **`supplementary_data`** contains supplementary data files referenced in the manuscript


**`scripts`** - scripts for data analysis

   - **`Clustifyr_classification.Rmd`** - Rmarkdown for cell sub type classification of murine clusters identified in `sc-DC-mmus.Rmd`

   - **`Geo_bulk_DeSEq2`** - Differential gene expression data of count-matrix in `res\GSE134869_SP7_umi_counts.tsv.gz'` `(1)`

   - **`Helft-anlysis`** - Differential gene expression analysis of cell sub types identified by Helft et al. `(2)` used for analysis in `Clustifyr_classification.Rmd`

   - **`sc-DC-mmus`** - Rmarkdown for analysis of host data

   - **`SC-Toxo`** - Rmarkdown for analysis of *Toxoplasma gondii* data

   - **`t.gondii_Crispr_analysis`** -Rmarkdown for comparative analysis using published data on in-vitro CRISPR-screen `(3)`

   - **`t.gondii_dc_function`** - Rscript containing functions used in Markdowns above 

## References

`(1) Pandey, S., Gruenbaum, A., Kanashova, T., Mertins, P., Cluzel, P., and Chevrier, N. (2020). Pairwise Stimulations of Pathogen-Sensing Pathways Predict Immune Responses to Multi-adjuvant Combinations. Cell Systems 11, 495-508.e10. 10.1016/j.cels.2020.10.001.`

`(2) Helft, J., Böttcher, J., Chakravarty, P., Zelenay, S., Huotari, J., Schraml, B.U., Goubau, D., and Reis e Sousa, C. (2015). GM-CSF Mouse Bone Marrow Cultures Comprise a Heterogeneous Population of CD11c+MHCII+ Macrophages and Dendritic Cells. Immunity 42, 1197–1211. 10.1016/j.immuni.2015.05.018.`

`(3) Wang, Y., Sangaré, L.O., Paredes-Santos, T.C., Hassan, M.A., Krishnamurthy, S., Furuta, A.M., Markus, B.M., Lourido, S., and Saeij, J.P.J. (2020). Genome-wide screens identify Toxoplasma gondii determinants of parasite fitness in IFNγ-activated murine macrophages. Nat Commun 11, 5258. 10.1038/s41467-020-18991-8.`


