library("metaseqR2")
current_directory <- getwd()
path1 <- file.path(current_directory)
#geneData <- loadAnnotation("mm10","ucsc","gene")
myTargets <- readTargets(file.path(current_directory, "targets_single.txt"))
###########################################################################################
result <- metaseqr2( sampleList = file.path(current_directory, "targets_single.txt"),
                     excludeList = NULL,
                     fileType = c("bam"),
                     contrast=c("KO_vs_WT"),
                     #  counts=mm9GeneCounts,
                     #  sampleList=sampleListMm9,
                     #  contrast=c("adult_8_weeks_vs_e14.5"),
                     #  libsizeList=libsizeListMm9,
                     refdb = "ucsc",
                     org = c("mm10"),
                     #  annotation="embedded",
                     #  embedCols=list(
                     #    idCol=4,
                     #    gcCol=5,
                     #    nameCol=8,
                     #    btCol=7
                     #  ),
                     countType="exon",
                     normalization="edger",
                     statistics = c("deseq", "edger", "nbpseq"), 
                     metaP = c("simes", "bonferroni", "fisher", "harmonic", "dperm.minp", "dperm.max", "dperm.weight", "fperm", "whitlock", "minp", "maxp", "weight", "pandora"), 
                     
                     pcut=0.05,
                     qcPlots=c(
                       "mds","filtered","correl","pairwise","boxplot","gcbias",
                       "lengthbias","meandiff","meanvar","deheatmap","volcano",
                       "mastat"
                     ),
                     figFormat=c("png","pdf"),
                     exportWhat=c("annotation","p_value","adj_p_value","fold_change","meta_p_value"),
                     exportScale=c("natural","log2"),
                     exportValues="normalized",
                     exportStats=c("mean","sd","cv"),
                     exportWhere=file.path(current_directory, "diff_exp"),
                     restrictCores=0.1,
                     geneFilters=list(
                       length=list(
                         length=500
                       ),
                       avgReads=list(
                         averagePerBp=100,
                         quantile=0.25
                       ),
                       
                       
                       expression=list(
                         median=TRUE,
                         mean=FALSE,
                         quantile=NA,
                         known=NA,
                         custom=NA
                       ),
                       biotype=getDefaults("biotypeFilter","mm10")
                     ),
                     outList=TRUE
)
###########################################################################################
