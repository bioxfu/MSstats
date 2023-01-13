'Usage:
  msstats.R --input=<DIR> --output=<DIR> --anno=<FILE> [--log]
  msstats.R (-h | --help)
  msstats.R --version

Options:
  --input=<DIR>    MaxQuant Input Folder 
  --output=<DIR>   MSstats Output Folder
  --anno=<FILE>    Experiment Designe Table
  --log            Keep Log Files
  -h --help        Show this screen.
  --version        Show version.
' -> doc

arguments <- docopt::docopt(doc, version = paste('MSstats', packageVersion('MSstats'), '\n'))
# print(arguments)

suppressPackageStartupMessages(library(MSstats))
input_protein_group <- paste0(arguments$input, '/proteinGroups.txt')
input_evidence <- paste0(arguments$input, '/evidence.txt')
input_annot <- arguments$anno
out_dir <-  arguments$output
out_dir <- paste0(out_dir, '/')
keep_logfile <- arguments$log

if (!dir.exists(out_dir)) dir.create(out_dir)

## First, get protein ID information
proteinGroups <- read.table(input_protein_group, sep = '\t', header = TRUE)

## Read in MaxQuant file: evidence.txt
infile <- read.table(input_evidence, sep = "\t", header = TRUE)

## Read in annotation including condition and biological replicates: annotation.csv
annot <- read.csv(input_annot, header = TRUE)

## 
quant <- MaxQtoMSstatsFormat(evidence = infile, annotation = annot, proteinGroups = proteinGroups, removeProtein_with1Peptide = TRUE)
maxquant.proposed <- dataProcess(quant, normalization = 'equalizeMedians',
                                 summaryMethod = "TMP",
                                 censoredInt = "NA", ## !! important for MaxQuant MBimpute=TRUE,
                                 maxQuantileforCensored = 0.999)
write.csv(quant, paste0(out_dir, '/MaxQtoMSstatsFormat.csv'), quote = FALSE, row.names = FALSE)
write.csv(maxquant.proposed$FeatureLevelData, paste0(out_dir, '/FeatureLevelData.csv'), quote = FALSE, row.names = FALSE)
write.csv(maxquant.proposed$ProteinLevelData, paste0(out_dir, '/ProteinLevelData.csv'), quote = FALSE, row.names = FALSE)

## 
dataProcessPlots(data = maxquant.proposed, type = "QCplot", ylimUp = 35,
                 width = 5, height = 5, which.Protein = 'allonly', address = out_dir)

dataProcessPlots(data = maxquant.proposed, type = "Profileplot", ylimUp = 35, 
                 featureName = "Peptide", width = 5, height = 5, address = out_dir)

dataProcessPlots(data = maxquant.proposed, type = "Conditionplot", 
                 width = 5, height = 5, address = out_dir)


# Comparing conditions with groupComparison
group_level <- levels(maxquant.proposed$ProteinLevelData$GROUP)
cbn <- t(combn(1:length(group_level), 2))[, 2:1]
comparison <- matrix(0, nrow = nrow(cbn), ncol = length(group_level))
rownames(comparison) <- paste(group_level[cbn[, 1]], group_level[cbn[, 2]], sep = '-')
colnames(comparison) <- group_level
for (i in 1:nrow(cbn)) {
  comparison[i, cbn[i, 1]] <- 1
  comparison[i, cbn[i, 2]] <- -1
}
maxquant.comparisons <- groupComparison(contrast.matrix = comparison, data = maxquant.proposed)
write.csv(comparison, paste0(out_dir, '/ComparisonMatrix.csv'), quote = FALSE)

# normal quantile-quantile plots
modelBasedQCPlots(data = maxquant.comparisons, type = "QQPlots", 
                  width = 5, height = 5, address = out_dir)
# residual plots
modelBasedQCPlots(data = maxquant.comparisons, type = "ResidualPlots", 
                  width = 5, height = 5, address = out_dir)

##
groupComparisonPlots(data = maxquant.comparisons$ComparisonResult, type = 'VolcanoPlot', 
                     width = 10, height = 10, address = out_dir)

pdf(paste0(out_dir, '/Heatmap_ColorKey.pdf'), width = 10, height = 4)
groupComparisonPlots(data = maxquant.comparisons$ComparisonResult, type = 'Heatmap', 
                     width = 10, height = 10, address = out_dir, colorkey = TRUE)
dev.off()

groupComparisonPlots(data = maxquant.comparisons$ComparisonResult, type = "ComparisonPlot", 
                     width = 5, height = 5, address = out_dir)

write.csv(maxquant.comparisons$ComparisonResult, paste0(out_dir, '/ComparisonResult.csv'), quote = FALSE, row.names = FALSE)

# Minimal number of biological replicates per condition
result.sample <- designSampleSize(data = maxquant.comparisons$FittedModel, 
                                  numSample = TRUE, desiredFC = c(1.25, 3), 
                                  FDR = 0.05, power = 0.8)
pdf(paste0(out_dir, '/SampleSize.pdf'))
designSampleSizePlots(result.sample)
dev.off()

if (!keep_logfile) file.remove(dir('.', 'MSstats_.*log'))

