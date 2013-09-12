##########################################################
# SETTINGS #
###########################################################
workDir <- "/Data/jyeung/projects/alternative_splicing/"
processedDataDir <- file.path("/Data/jyeung/projects/alternative_splicing/output/tables")
outputDir <- file.path(workDir, "output/optdis_output")

scriptDir <- "/Data/kwang/bdvtools/TestFramework"
snDir <- "/Data/kwang/bdvtools/Subnetwork"

#. Set parameters and load libraries .#
options(width = 100, digits=8)
source(file.path("/Data/kwang/bdvtools/ArrayProcessing", "visualizeBatch.R"))
source(file.path(scriptDir, "runAnalysis.R"))
source(file.path(scriptDir, "utilities.R"))
###########################################################

file.input <- file.path(processedDataDir, "gene_expression_data_forOptDis_samplenames.txt")

dat <- read.table(file.input, sep="\t", row.names=1, header=TRUE, stringsAsFactors=F)

Train_expr <- dat
Test_expr <- dat


#Note: execute the runExternalValidation_SN() command as given below and cancel the execution immediately. This will create the directory structure as given below
#put the training and test expression matrix in the following path: /Data/jyeung/projects/alternative_splicing/output/optdis_output/activity_gene/MarkerDiscovery/takeda/subtype/
#path: outputDir/activity_gene/analysis/batch/endpoint/
#rename Training and Test expression matrix as Training.txt and Test.txt respectively
#after placing the required files in the path given run runExternalValidation_SN() command to its full execution

source(file.path(scriptDir, "runAnalysis.R"))
runExternalValidation_SN(trainExpr=Train_expr, testExpr=Test_expr, analysis="MarkerDiscovery", batch="takeda", endpoint="subtype", 
	network="HPRD", optSN="OptDis", optActivity="Single", scale="Brief", platform="U133Plus2", sameplatform=T, ftrSelection=F, supervisedLearning="Classification")
	
plotPerformance_CustomSN(analysis="MarkerDiscovery", batch="takeda", endpoint="subtype", optSN="OptDis", optActivity="Single")


