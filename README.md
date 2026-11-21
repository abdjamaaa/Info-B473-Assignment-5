# Assignment 5: Advanced R

Abdul Djama  
INFO-B473  
November 21, 2026

## What I Did

This assignment was done in R to compare the expression of genes between the tumor and the normal. I began by sorting out the data and pairing the samples with the appropriate labels. Next, I separated the tumor and normal samples, took the mean of the gene expression and determined the fold changes. Those results I used for making graphs to show the differences between the 2 groups.

## Files Used

- Gene_Expression_Data.xlsx
- Gene_Information.csv
- Sample_Information.tsv
- Info-B473 Assignment5.R

## How to Run the Code

The R script and the three data files must be placed in the same directory to run the code. Open the script in RStudio and change the working directory to that folder. After the code has executed it will generate a folder called results, which will include the graphs and the data tables.

These packages will be required: readxl, ggplot2, dplyr, tidyr and pheatmap.

If the packages aren't installed, use this command:

`install.packages(c("readxl", "ggplot2", "dplyr", "tidyr", "pheatmap"))`

## Results

The dataset had 48,803 probes and 18 samples. There were nine samples of tumors and nine samples of normal tissue. Once I calculated the fold changes, I got 7 probes with the absolute log2 fold change > 5.

Of the 7 probes, 4 were expressed at higher level in tumor tissue and 3 were expressed at higher level in normal tissue. This yielded 57.1% higher in tumour and 42.9% higher in normal.

There were 3 probes that were differentially expressed on chromosome 2.

I plotted 3 bar graphs on the chromosome count, the difference between the chromosome count of the tumor and normal samples and the percentage of higher expression in each group. I also produced a heatmap and a clustermap to see how the gene expression is different in different samples.

On the heatmap and clustermap, I plotted the 100 probes that most varied to make the patterns easier to see.

## Fold Change Calculation

The formula given in the assignment was: log2((Tumor - Normal) / Normal). But this formula does not apply when the expression of the tumor is lower than or equal to normal expression.

With this in mind, I also calculated the ratio of gene expression as normal over tumor and tumor over normal, and performed these calculations using `log2(Tumor / Normal)` and `log2(Normal / Tumor)`. Results were based on the fold change cutoff, not on statistical significance. For the 7 selected results, these are probes, and not necessarily different genes.

## Saved Results

The code saves the following files in the results folder:

- 01_expression_renamed.csv
- 02_tumor_data.csv
- 03_normal_data.csv
- 04_fold_changes.csv
- 05_selected_degs.csv
- 06_chromosome_counts.png
- 07_chromosome_groups.png
- 08_expression_percentages.png
- 09_heatmap.png
- 10_clustermap.png
- 11_summary.txt
