<img width="600" alt="Screenshot 2022-02-12 at 18 15 43" src="https://user-images.githubusercontent.com/63663530/153721166-e192850b-f449-4530-8166-dc7faa79790a.png">

Developed during my academic internship at Instituto de Investigación Sanitaria fundación Jiménez Díaz under the tutelage of Dr Pablo Mínguez Paniagua and Dra Gloria Álvarez-Llamas

# What is EnrichProt? 


**EnrichProt** is an integral tool that allows researchers to extract the underlying biological information among the vast amount of proteomics data. It is presented as a utility that helps users evaluate and select certain parameters in their analysis by studying their influence on biological information.

The two main parameters that can be evaluated are:
1. The fold change
2. The number of peptides with which a protein has been identified

### Fold change

The first parameter EnrichProt can evaluated is the fold change, used as a threshold to detect significantly altered proteins in a condition of study compared to a control. It gives information about the degree of change between the two conditions or biological classes under study. However, determining the optimum value is often complicated as it will depend on each experiment. Too strict or too relaxed fold change values would lead to misleading outcomes. 

### Number of peptides 

The second parameter EnrichProt can evaluate is the most appropriate number of peptides that identifies a protein to be considered. The study of the number of peptides with which a protein has been identified has never been reported in any tool or method. In MS-based proteomics analysis, proteins can be detected with a single peptide, which complicates the analysis due to the lack of knowledge of whether it is a real hit or just a false positive as a result of a contamination during sample handling.

### EnrichProt main objective

In order to elucidate which criteria are the most appropriate for the analysis, EnrichProt evaluates several fold changes and number of peptides simultaneously, extracting the relevant biological information from the vast amount of detected proteins for each parameter selection. This analysis allows us to analyze which criteria reports a signal with less noise, due to fewer false positives, and which criteria offers a stronger signal due to fewer false negatives. In addition to provide detailed analysis for each of the parameters under study, EnrichProt provides the researcher an overview of the obtained results and of the analysis performance.

# EnrichProt framework

Three types of analysis can be performed using EnrichProt pipeline: 

1) **Differential abundance analysis**
2) **Functional profiling using GSEA**
3) **Functional profiling using ORA**

The study of the influence of the number of peptides and the fold change can be performed automatically selecting several values for each one. The general framework of EnrichProt pipeline is shown below:



![image](https://user-images.githubusercontent.com/63663530/153719183-356ecb98-70fb-442a-b3b6-9c15f8e93e11.png)


## <img width="52" alt="Screenshot 2022-02-12 at 22 24 21" src="https://user-images.githubusercontent.com/63663530/153728996-7908aa06-6e11-4bdb-9bdb-6c56dc96b4f8.png"> Differential abundance analysis




In order to detect changing proteins in the condition of interest, a **Welch's two-sample t-test** followed by a multiple testing correction is performed by **FDR**. A simple **annotation** of proteins using Gene Ontology is also performed.

### Input:

+ EnrichProt receives a file containing identified UniProt IDs (1st column) and the different values for each feature in each of the samples under study. 
+ To study the influence of the number of peptides, a column containing this value is required. 
+ The two groups under study with their column positions are needed together with the threshold values of the parameters to evaluate (p-value, fold change and number of peptides). 
+ The user can choose between calculating the fold change or the log2(fold change).  

### Output:

+ EnrichProt returns the original input matrix where calculated statistical parameters (t-statistic, fold change, p-value and adjusted p-value by FDR) and protein annotations are added. 
+ Multiple graphics are generated and the same results are also provided for detected proteins with each selected fold change if several values have been chosen.  

<img width="953" alt="Screenshot 2022-02-12 at 20 07 04" src="https://user-images.githubusercontent.com/63663530/153724831-ff3f2619-1cee-44c4-a1db-487dc9fb354d.png">



## <img width="53" alt="Screenshot 2022-02-12 at 22 23 56" src="https://user-images.githubusercontent.com/63663530/153729035-143eba73-b2ab-4647-9eab-a95e8d7ecd96.png"> Functional profiling 

**Functional profiling methods** allow the interpretation of proteomics data in terms of biological functions, pathways, and protein interactions that are deregulated at a particular cellular condition. This kind of analysis detects statistically over-represented modules of genes or proteins with a functional entity in the complete set of identified elements in case-control experiments or similar approaches. 

Two popular methods are **Over Representation Analysis (ORA)** and **Gene Set Enrichment Analysis (GSEA)**, both using pre-defined sets of genes that share common biological functions or cell locations as modules to test.


## <img width="26" alt="Screenshot 2022-02-12 at 22 50 11" src="https://user-images.githubusercontent.com/63663530/153729910-c3d2dc92-5758-4318-bc25-b5c7d7055a94.png"> Functional profiling using ORA

ORA requires a **threshold** that select significantly deregulated proteins compared to the rest of the list of identified proteins (background). The choice of this threshold is often a crucial decision that depends on the number of observations being tested.

### Input:

+ EnrichProt receives a file containing identified UniProt IDs (1st column) and a column with the parameter chosen to represent a phenotype variation. 
+ In a case-control experiment, this parameter would be the statistics or the fold change, which can be calculated by EnrichProt in the differential abundance analysis or provided by the user. EnrichProt ranks the dataset according to the chosen parameter.
+ In ORA the **selection of a threshold** to define significantly changing proteins (e.g., fold change) is required. Multiple values of the selected parameter can be chosen in order to perform an analysis on the influence of each value in the obtained information. 
+ The adjusted p-value can be selected to retain, from all enriched categories, those below these values. 
+ Gene sets can be filtered out specifying the minimal and maximal size (number of annotated proteins) to be considered, being 10 and 500 the default ones. 
+ Ontology sets that can be used: GO, KEGG, DisGeNET, Molecular Signatures Database (MSigDB) and The Human Protein Atlas (HPA). 
+ In the case of selecting GO, EnrichProt provides all enriched GO terms and an additional list where redundant terms have been removed, performing the Relevance semantic similarity method in *REVIGO* package and using a small similarity threshold (C=0.5). 

### Output:

+ The output consists of significantly enriched functional categories under the selected p-value for each peptide, database and threshold.
+ Several informative figures are generated. 

![image](https://user-images.githubusercontent.com/63663530/153730169-c8983633-4eed-4c97-8ac7-a505c1024b4c.png)


## <img width="26" alt="Screenshot 2022-02-12 at 23 02 37" src="https://user-images.githubusercontent.com/63663530/153729915-29ee71e3-7a51-4022-8425-f2220e652bb9.png"> Functional profiling using GSEA

GSEA attempts to measure the coordinated behavior of genes/proteins from a **ranked list** based on a parameter linked to the phenotype we want to describe. GSEA **does not require a threshold** to define significantly abundant proteins like ORA, and uses the distribution of pre-defined gene sets within the ranked list to determine whether that category is up or down-regulated.

### Input:

+ EnrichProt receives a file containing identified UniProt IDs (1st column) and a column with the parameter chosen to represent a phenotype variation. 
+ In a case-control experiment, this parameter would be the statistics or the fold change, which can be calculated by EnrichProt in the differential abundance analysis or provided by the user. EnrichProt ranks the dataset according to the chosen parameter.
+ The adjusted p-value can be selected to retain, from all enriched categories, those below these values. 
+ Gene sets can be filtered out specifying the minimal and maximal size (number of annotated proteins) to be considered, being 10 and 500 the default ones. 
+ Ontology sets that can be used: GO, KEGG, DisGeNET, Molecular Signatures Database (MSigDB) and The Human Protein Atlas (HPA). 
+ In the case of selecting GO, EnrichProt provides all enriched GO terms and an additional list where redundant terms have been removed, performing the Relevance semantic similarity method in *REVIGO* package and using a small similarity threshold (C=0.5). 

### Output:
+ The output consists of significantly enriched functional categories under the selected adjusted p-value for each selected peptide and ontology. 
+ Several informative figures are generated. 

![image](https://user-images.githubusercontent.com/63663530/153730090-03b6b4f7-bbec-4f59-acd1-490fea07dd28.png)








