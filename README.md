# Parental effects in *Neurospora crassa*
Data and scripts for the manuscript "Parental effects in a filamentous fungi: phenotype, fitness and mechanism" by Mariana Villalba de la Pena, Pauliina A.M Summanen, Neda N. Moghadam and Ilkka Kronholm

## **STATISTICAL ANALYSIS**
#### The file stat_analysis.R contains R scripts for the analysis of:
- Initial growth of the two generations, separated by experiments and all the experiments                       together.
- Growth rate.
- Initial growth rate of the mutant strains.
- Number of conidia produced in 1.5% and 0.015% sucrose environment.
- Conidia viability.
- Conidia size.
- Initial growth in alternative carbon sources.
- Protein and carbonhydrate content (i.e glucose, glycogen)

#### The data files:
- growth_measure_final_plating.csv : mycelium growth database of both generations for all                  experiments. 
- sorbose_viability.csv : number of colonies growing in sorbose plates as a viability proxy.
- tgen_mutants : mycelium growth of mutant strains.
- number_conidia.csv : number of conidia produce in each of the sucrose environment.
- protein_assay020822.csv : protein amount in conidia in each of the sucrose environment.
- glucose_assay141022.csv : glucose content in conidia in each of the sucrose environment.
- glycogen_assay120822.csv : glycogen content in conidia in each of the sucrose environment.

## **Competition experiments**
#### Script files for analysis are in folder compexp
- compexp.R : contains R scripts to analyze competition experiments
- HRMqPCR.R : contains R scripts required for handling of HRM-PCR data
#### Data files
- compexp1.RData : contains processed data for the first set of competition experiments
- compexp2.RData : contains processed data for the second set of competition experiments
- 190812_hrm.csv : raw data for plate 1 set 1
- 190812_samples.csv : sample information for plate 1 set 1
- 190812_2_hrm.csv : raw data for plate 2 set 1
- 190812_2_samples.csv : sample information for plate 2 set 1
- 190812_3_hrm.csv : raw data for plate 3 set 1
- 190812_3_samples.csv : sample information for plate 3 set 1
- 190917_PS_hrm.csv : raw data for plate 1 set 2
- 190917_PS_samples.csv : sample information for plate 1 set 2

## **RNA-seq**
#### Hisat bash script.
#### In the RNA_seq.R script is the following code for RNA-seq data analysis:
- Remove unwanted variation using RUVseq bioconductor package (R environment version 4.0.2) using              positive and negative ERCC spike-in controls.
- Differential gene expression analysis using edgeR and DESEq2.
- KEGG pathways enrichment analysis (ORA and GSEA) using cluster profiler.
- Count number of genes belogning to enrich GSEA enrich pathways.
#### The data files:
- all_hista_counts_mod.txt : counts, hisat output.
- genes_path.csv : number of genes in each pathway.
