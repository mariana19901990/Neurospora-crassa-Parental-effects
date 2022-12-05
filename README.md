# Neurospora-crassa-Parental-effects
Data scripts and script of manuscript "Parental effects in a filamentous fungi: phenotype, fitness and mechanism"
# N.crassa-parental-effects
Data for manuscript "Parental effects in a filamentous fungus: phenotype, fitness and mechanism" by Mariana Villalba de la Pena, Pauliina A.M Summanen, Neda N. Moghadam and Ilkka Kronholm

## **STAT ANALYSIS**
#### In the stat_analysis.R script code for the analysis of:
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
