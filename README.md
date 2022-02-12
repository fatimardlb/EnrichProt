# <img width="698" alt="Screenshot 2022-02-12 at 18 15 43" src="https://user-images.githubusercontent.com/63663530/153721166-e192850b-f449-4530-8166-dc7faa79790a.png">

An integral pipeline to extract relevant biological information among the vast amount of proteomics data.

## What is EnrichProt?

**EnrichProt** is an integral tool that allows researchers to extract the underlying biological information among the vast amount of proteomics data. It is presented as a utility that helps users evaluate and select certain parameters in their analysis by studying their influence on biological information.

The two main parameters that can be evaluated are:
1. The fold change
2. The number of peptides with which a protein has been identified

### Fold change

The first paraemter EnrichProt can evaluated is the fold change, used as a threshold to detect significantly altered proteins in a condition of study compared to a control. It gives information about the degree of change between the two conditions or biological classes under study. However, determining the optimum value is often complicated as it will depend on each experiment. Too strict or too relaxed fold change values would lead to misleading outcomes. 

### Number of peptides 

The second parameter EnrichProt can evaluate is the most appropriate number of peptides that identifies a protein to be considered. The study of the number of peptides with which a protein has been identified has never been reported in any tool or method. In MS-based proteomics analysis, proteins can be detected with a single peptide, which complicates the analysis due to the lack of knowledge of whether it is a real hit or just a false positive as a result of a contamination during sample handling.

### EnrichProt main objective

In order to elucidate which criteria are the most appropriate for the analysis, EnrichProt evaluates several fold changes and number of peptides simultaneously, extracting the relevant biological information from the vast amount of detected proteins for each parameter selection. This analysis allows us to analyze which criteria reports a signal with less noise, due to fewer false positives, and which criteria offers a stronger signal due to fewer false negatives. In addition to provide detailed analysis for each of the parameters under study, EnrichProt provides the researcher an overview of the obtained results and of the analysis performance.

## EnrichProt framework

![image](https://user-images.githubusercontent.com/63663530/153719183-356ecb98-70fb-442a-b3b6-9c15f8e93e11.png)





