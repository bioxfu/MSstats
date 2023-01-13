# MSstats

## 1. Build Docker image
cd docker
./build.sh

## 2. Run
./run.sh

## 3. Output
### Quality control and normalization effects
**QC plot** (*QCPlot.pdf*) visualizes potential systematic biases between mass spectrometry runs. 
It can be used to assess the effects of the normalization step. After constant normalization, the median intensities of reference transitions across all proteins should be equal between runs. After quantile normalization, the distribution of reference intensities across all proteins should be equal between runs. 

**Profile plot** helps identify potential sources of variation for each protein. 
Profile plots with summarization present the effects of the summarization step by showing all individual measurements of a protein and their summarized intensity. 
The first plot (*ProfilePlot.pdf*) shows individual measurement for each peptide (peptide for DDA, transition for SRM or DIA) across runs, grouped per condition. Each peptide has a different color/type layout. Disconnected lines show that there are missing value (NA). 
The second plot (*ProfilePlot_wSummarization.pdf*) shows run-level summarized data per protein. The same peptides (or transition) in the first plot are presented in grey, with the summarized values (by TMP, in this example) overlaid in red. 

**Condition plot** (*ConditionPlot.pdf*) visualizes potential systematic differences in protein intensities between conditions. Dots indicate the mean of log2 intensities for each condition. With the option interval='CI'(default), error bars indicate the confidence interval with 0.95 significant level for each condition. With the option interval='SD', error bars indicate the standard deviation among all feature intensities for each condition. 

### Verifying the assumption of the model
**Normal quantile-quantile plot** (*QQPlot.pdf*) illustrates that such deviations from constant variance can be mistaken for deviations from Normality. Only large deviations of transition intensities from the straight line are problematic.

**Residual plot** (*ResidualPlot.pdf*) shows variance of the residuals that is associated with the mean feature intensity. Any specific pattern, such as increasing or decreasing by predicted abundance, is problematic.

### Visualization of differentially abundant proteins
**Volcano plots** (*VolcanoPlot.pdf*) visualize the outcome of one comparison between conditions for all the proteins.
The y-axis displays the FDR-adjusted p-values on the negative log10 scale, representing statistical significance. The horizontal dashed line shows the FDR cutoff. 
The points above the FDR cutoff line are statistically significant proteins that are differentially abundant across conditions. These points are colored in red and blue for upregulated and downregulated proteins, respectively. 
The x-axis is the model-based estimate of fold change on log scale, and represents practical significance. 
If the fold change cutoff is specified, the points above the horizontal cutoff line but within the vertical cutoff line will be considered as not differentially abundant (and will be colored in black). 

**Heatmaps** (*Heatmap.pdf*, *Heatmap_ColorKey.pdf*) illustrate the patterns of up- and down-regulation of proteins in several comparisons. 
Columns in the heatmaps are comparison of conditions assigned in contrast.matrix, and rows are proteins. 
The heatmaps display signed FDR-adjusted p-values of the tests, colored in red/blue for significantly up- /down-regulated proteins. Brighter colors indicate stronger evidence in favor of differential abundance. Black color represents proteins that are not significantly differentially abundant. 
The rows and columns of the heatmaps can be ordered with the option clustering, which performs hierarchical clustering with the Ward method (minimum variance). 

**Comparison plots** (*ComparisonPlot.pdf*) illustrate model-based estimates of log-fold changes, and the associated uncertainty, in several comparisons of conditions for one protein. 
X-axis is the comparison of interest. Y-axis is the log fold change. The dots are the model-based estimates of log-fold change, and the error bars are the model-based 95% confidence intervals (the option sig can be used to change the significance level of significance). 
For simplicity, the confidence intervals are adjusted for multiple comparisons within protein only, using the Bonferroni approach. For proteins with N comparisons, the individual confidence intervals are at the level of 1-sig/N. 

### Sample size calculation for a future experiment
**Sample size calculation** (*SampleSize.pdf*) The calculated relationship between the number of biological replicates per condition (numSample), average statistical power across all the proteins (power), minimal fold change that we would like to detect (desiredFC), and the False Discovery Rate (FDR) can be visualized using the function designSampleSizePlots. The function takes as input the output of designSampleSize. 


### CSV files
- The output of MaxQtoMSstatsFormat(): MaxQtoMSstatsFormat.csv	
- Output of dataProcess(): FeatureLevelData.csv and ProteinLevelData.csv
- Contrast matrix used by groupComparison(): ComparisonMatrix.csv
- Output of the groupComparison(): ComparisonResult.csv
