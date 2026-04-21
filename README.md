# 1. [findseq.py](./findseq/findseq.py)
## System requirement
Runs on Python3. Tested on Python 3.8.10. No non-standard hardware required.
## Installation guide
No installation required. Grant execution permission by `chmod +x findseq.py`
## Demo
Find demo input under findseq/T2T-CHM13v2.0.fna (hs1 genome) and run via `python3 findseq.py T2T-CHM13v2.0.fna CpG.bed CG`. Output should be locations of all occurrences of CpG (CG) on hs1 genome in BED format. <10 min. runtime is expected.
## Instructions for use
Find dinucleotide (or sequence motif of any length) via `python3 findseq.py <Genome.fasta> <Out.bed> <Motif>`.

# 2. [bed_2_tbtools.py](./bed_2_tbtools/bed_2_tbtools.py)
## System requirement
Runs on Python3. Tested on Python 3.12.8. No non-standard hardware required.
## Installation guide
No installation required. Grant execution permission by `chmod +x bed_2_tbtools.py`
## Demo
Find demo input under bed_2_tbtoos/CG.bed (All CpG sites of hs1 chromosome 1). Run script by `python3 bed_2_tbtools.py`. Enter input (`CG.bed`) and output (`CG.BINstat.tab.xls`) filenames when prompted. A file summarizing the length of region occupied by the BED file in every 10,000 bp window of the genome is expected. <10 min. runtime is expected.
## Instructions for use
Run script while inputting your BED file at prompt. The output is compatible as input for visualization as a circos plot track in TBtools ([Chen et al., 2023](https://www.cell.com/molecular-plant/fulltext/S1674-2052(23)00281-2))

# 3. [intro.R](./intro/intro.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Execute directly on R by blocks specified in the code, preferrably on Rstudio, by importing each input file uploaded under [intro](./intro) to R session as df or the variable name used in each block via
```
df=read.table("<Filename>")
```
e.g. 
```
df=read.table("Fig1b.tsv")
```
and prompting the block. A visualization of figure in the format of Fig. 1b-g of paper ([Lee et al., 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected. Less than 5 minute run time expected for each block.
## Instructions for use
Import your data via
```
df=read.table("<Filename>")
```
replacing `<Filename>` with your data before executing a block. The data must be a tab-delimited file with column names specified in each block.

# 4. [plot_profile.R](./plot_profile/plot_profile.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Declare functions `preprocess_data()` and `create_tss_plot()` from script. Pipe the demo input found under plot_profile/MP_TSS_hg38.tsv through functions: I. `preprocess_data()`, II. `create_tss_plot()` and III. print output
e.g.
```
processed_data <- preprocess_data("MP_TSS_hg38.tsv")
tss_plot <- create_tss_plot(processed_data)
print(tss_plot)
```
A line plot of an average MP profile at TSS as in format of Fig. 2 of study ([Lee et al., 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1.full)) is expected. <5 min. runtime is expected.
## Instructions for use
Prepare a tab-delimited file with columns indicating information on I. Vicinal TSS ID, II. Distance from TSS, III. MP and IV. Strand (+/-) for every row corresponding to a CpG. Input the file through the aforementioned pipeline:
```
processed_data <- preprocess_data("<Your_file.tsv>")
tss_plot <- create_tss_plot(processed_data)
print(tss_plot)
```

# 5. [compute_tvd.R](./compute_tvd/compute_tvd.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Declare function `compute_tvd_pvalues()` from script. Find inputs compute_tvd/MP_wg_hg38.tsv (MP of all CpGs on hg38 genome) and compute_tvd/MP_TSS_hg38.tsv (MP and distance from TSS for all TSS-vicinal (±10,000 bp) CpGs on hg38). Run function via
```
compute_tvd_pvalues("MP_wg_hg38.tsv", "MP_TSS_hg38_tsv")
```
A data frame comprising information on p-value at every distance point is expected. 30 hour runtime is expected.
## Instructions for use
Prepare two tab-delimited files: 1. Population (`wg_file`; whole-genome) and 2. sample (`element_file`; region of interest e.g. TSS vicinity). The population file must comprise rows corresponding to every CpG on genome and contain a column with MP information (column 4 by default). The sample file must comprise rows corresponding to CpG inside the region of interest and contain a column with distance, MP and strand information (columns 2, 3 and 4, respectively, by default). Run function via
```
compute_tvd_pvalues("<population_file>", "<sample_file>")
```

# 6. [visualize_tvd.R](./visualize_tvd/visualize_tvd.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Start R session and declare function `plot_pvalue_heatmap()` from script. Import demo input found at [visualize_tvd/pvalues_TSS_hg38.tsv](./visualize_tvd/pvalues_TSS_hg38.tsv) via
```
df=read.table("pvalues_TSS_hg38.tsv")
```
Run function with input via
```
plot_pvalue_heatmap(df, "Title")
```
A horizontally linear heatmap indicating p-value at each distance through color intensity as in Fig. 2 of ([Lee et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected. <5 min. runtime expected.
## Instructions for use
Prepare a data frame format file with two columns denoting distance from a site of interest and p-value at the distance, preferably an output from [compute_tvd.R](./compute_tvd/compute_tvd.R). Run function via
```
plot_pvalue_heatmap(<p_value_heatmap>, "<Your_title>")
```

# 7. [scatterplot_gene.R](./scatterplot_gene/scatterplot_gene.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Start R session and declare function `plot_gene_mp()` from script. Import demo input found under [scatterplot_gene/scatterplot_sample.tsv](.scatterplot_gene/scatterplot_sample.tsv) via
```
df=read.table("scatterplot_sample.tsv")
```
Run function via
```
plot_gene_mp(df, "Class", c("Mammalia"))
```
A scatterplot displaying MP of all CpG of mammals at TSS of gene ACTB, as an example, as in format of Extended Data Fig. 7c in study ([Lee et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected. <5 min. runtime is expected.
## Instructions for use
Prepare a data frame format file with columns Distance, MP and Class indicating distance and MP of CpG and which class it belongs to, respectively. Run function indicating input file, the column indicating class information and which class you want to visualize i.e.
```
plot_gene_mp(<Data_frame>, <Class_column_name>, <Class_of_interest>)
```
Class must be one or more of: `Mammalia`, `Aves`, `Reptilia`, `Amphibia`, `Sarcopterygii`, `Actinopterygii` and `Chodrichthyes`.



# 9. [plot_profile_promoter.R](./plot_profile_promoter/plot_profile_promoter.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Start R session and declare function `plotMedian()` from script. Feed demo input found under [plot_profile_promoter/medianMP_transcriptTSS_T2Thuman.tsv](./plot_profile_promoter/medianMP_transcriptTSS_T2Thuman.tsv) (Median MP profile at TSS of hs1) to function while specifying promoter start, end, desired line color and core promoter length by
```
plotMedian("medianMP_transcriptTSS_T2Thuman.tsv", -874, 948, "#000000", 170)
```
A median MP profile with the promoter and core promoter size indicated as in Fig. 5f,g and Fig. 6 of study ([Lee et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected.
## Instructions for use
A tab-delimited file summarizing median methylation probability (MP) values at every distance within ±10,000 bp of site of interest must be calculated beforehand and provided as input. Run function via
```
plotMedian("<Filename>", <Promoter_start>, <Promoter_end>, <Color(Hex)>, <Core_size>)
```
to visualize the median MP profile of species.

# 10. [plot_comparison.R](./plot_comparison/plot_comparison.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Start R session and preceed running whole script by importing demo input found under [plot_comparison/comparison_sample.tsv](./plot_comparison/comparison_sample.tsv) as `df` via
```
df=read.table("comparison_sample.tsv")
```
A scatterplot figure in the format of Fig. 5h-j of study ([Lee et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected.
## Instructions for use
Organize your calculated promoter lengths for individual species beforehand into a data frame format with columns `Class`, `Genome.size` and `Promoter.size` and preceed prompting of the script with input importing:
```
df=read.table("<Your_promoter_sizes.tsv>")
```

# 11. [plot_umap.R](./plot_umap/plot_umap.R)
## System requirement
Runs on R. Tested on Rstudio/2025.09.2+418 "Cucumberleaf Sunflower" Release powered by R version 4.4.1. No non-standard hardware required.
## Installation guide
No installation required.
## Demo
Start R session and declare function `plot_UMAP()`. Import demo input found under [plot_umap/umap_sample.tsv](./plot_umap/umap_sample.tsv) via
```
df=read.table("umap_sample.tsv")
```
and conduct visualization via
```
plot_UMAP(df, "Title", df$Class, class_color)
```
A scatterplot summarizing UMAP components in a two-dimensional space as in the format of Fig. 4 of study ([Lee et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.29.715150v1)) is expected.
## Instructions for use
Organize UMAP result into a data frame format with columns `UMAP1` and `UMAP2` for the top two UMAP components. Specify column name containing information color coding (e.g. phylogenetic class or tissue type). Specify a vector defining color coding scheme. Run via
```
plot_UMAP(<Data_frame>, <Title>, <Column_name>, <Color_coding_scheme>)
```
