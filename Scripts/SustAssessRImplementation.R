#Script for implementation of SustAssessR pipeline
#Author: Nick Matthews

######Set-up#####
library(stringr)
library(EnvStats)
library(magrittr)
library(optiRum)

#####Directory set-up#####
#set data directory
parameter.dir <- file.path("Data/Parameters")

#Set output directory
output.dir <- "Data/Results"
data.dir <- "Data"
dir.create(file.path(output.dir),showWarnings = FALSE)

#Need to load in necessary source code
source("Scripts/SustAssessRSourceCode.R")

#Load in modules
source("Scripts/ProcessModules.R")


#####Set initial parameters#####
#This section specifies the key parameters for the modelling
#Higher-level parameters
FU <- 1 #kg product
targets = c("Cadaverine","Putrescine")
perspectives <- c("Hierarchist","Egalitarian","Individualist")
sugarScen = c(BR="sucrose",US="glucose",USLigno="Mixed",IN="sucrose",FR="sucrose",DEStraw="glucose",DEOak="glucose")

#Parameters for dataframe
NitrogenSources <- c("AmmoniumSulfate","AmmoniumNitrate")
energySources <- c("Grid","Biomass")
WasteScenario<-c("Non-Integrated","Integrated")
AdipicAcidScens <- c("AdipicAcid_RoW","AdipicAcid_Reduct")
SebacicAcidScens <- c("SebacicAcid_BestCaseCN","SebacicAcid_WorstCaseCN")
geogs = c("BR","US","USLigno","FR")
DSPs <- c("80","95")
Switch <- c("0hours","6hours","12hours")
OrgSpeed <- c("Medium","Fast")
MaxYield <- c("50%","75%","90%")
times <- c(10:30)

#Mpnte-carlo runs
MCruns = 10000
#Specify number of multi-starts for sensitivity test
MultiStart = 1000

#Name of the run
overallName <- paste0(paste(targets,collapse=""))

#####Load in data#####
#Specify links
linkages <- list()
linkages[["linkagesNon-Integrated"]] <- read.csv(file.path("Data/Linkages","LinkagesNonIntegrated.csv"),col.names=c("type","from","to","flow"),stringsAsFactors = FALSE)
linkages[["linkagesIntegrated"]] <- read.csv(file.path("Data/Linkages","LinkagesIntegrated.csv"),col.names=c("type","from","to","flow"),stringsAsFactors = FALSE)
linkages[["linkagesNylonPolymerisation"]] <- read.csv(file.path("Data/Linkages","LinkagesNylonPolymerisation.csv"),col.names=c("type","from","to","flow"),stringsAsFactors = FALSE)


#Load in parameters
allParameters <- read.csv(file.path(parameter.dir,"ParameterSet.csv"),col.names=c("Par","Value","Min","Max","SD","Dist","Units","Notes","Type"),row.names=1,stringsAsFactors = FALSE)
sensTestPars <- read.csv(file.path(parameter.dir,"SensTestParameters.csv"),stringsAsFactors = FALSE)

#Set up the original parameters data.frame
discreteparameters <- data.frame(geogs=sample(geogs,MCruns,replace=TRUE),DSPs=sample(DSPs,MCruns,replace=TRUE),
                         Switch=sample(Switch,MCruns,replace=TRUE),OrgSpeed=sample(OrgSpeed,MCruns,replace=TRUE),
                         MaxYield=sample(MaxYield,MCruns,replace=TRUE),WasteScenario=sample(WasteScenario,MCruns,replace=TRUE),
                         times=sample(times,MCruns,replace=TRUE),energySource=sample(energySources,MCruns,replace=TRUE),
                         NitrogenSource=sample(NitrogenSources,MCruns,replace=TRUE),
                         SebacicAcidComp=sample(SebacicAcidScens,MCruns,replace=TRUE),AdipicAcidComp=sample(AdipicAcidScens,MCruns,replace=TRUE),
                         stringsAsFactors = FALSE)
vardists <- makeDists(allParameters,samples = MCruns)
parameters <- cbind(discreteparameters,vardists)
#Load in or create the ReCiPe characterisation factors
CFs <- recipeCFs("Data/ReCiPe")

#Load in background data for all scenarios
backgroundfiles <- dir(file.path("Data/Background"))[grepl(".csv",dir(file.path("Data/Background")))]
backgroundData <- list()
for(filename in backgroundfiles){
  temp <- read.csv(file.path("Data/Background",filename),header = FALSE,na.strings = c("","NA"),stringsAsFactors = FALSE)
  temp <- cleanGabiResults(temp)
  rownames(temp) <- temp$Indicator
  temp$Indicator <- NULL
  backgroundData[[gsub(".csv","",filename)]] <- temp
}

#Load in all the process flows
flowfiles <- dir(file.path("Data/ProcessedFlows"))
flowData <- list()
flowRefsDF <- read.csv(file.path("Data/Parameters","FlowRefs.csv"),row.names=1,stringsAsFactors = FALSE)
for(filename in flowfiles) {
  temp <- read.csv(file.path("Data/ProcessedFlows",filename),stringsAsFactors = FALSE)
  colnames(temp) <- c("AllDetails","Unit","Name","Input_output_parameter",gsub(" ","",flowRefsDF$Name))
  flowData[[gsub(".csv","",filename)]] <- temp
}

#Load in comparisons file
comparisons.dir <- file.path("Data/Comparisons")
comparisonName <- "NylonComparisons"
perspectives <- c("Egalitarian","Hierarchist","Individualist")
comparisonFiles <- list()
for(ii in 1:length(perspectives)){
  tempfile <- read.csv(file.path(comparisons.dir,paste0(comparisonName,"_",perspectives[ii],".csv")),header = FALSE,na.strings = c("","NA"),stringsAsFactors = FALSE)
  tempfile <- cleanGabiResults(tempfile)
  row.names(tempfile) <- tempfile$Indicator
  comparisonFiles[[paste0("Comparisons",perspectives[ii])]] <- tempfile
}


#####Run the analysis#####
#Calculate the results
results <- allResultsProcess(data.dir,linkages,parameters,flowData,backgroundData,CFs,FU=FU,targets=targets,perspective=perspectives,MCruns=nrow(parameters))

#Normalise the results
resultsNormalised <- normalisation(data.dir="Data/ReCiPe",results,perspectives=c("Egalitarian","Hierarchist","Individualist"),targets=c("Cadaverine","Putrescine"))

save(results,resultsNormalised,
     file=file.path(output.dir,paste0("NormalisedResults",overallName,".RData")))

#Can only do hierarchist for adipic acid sensitivity case of reduced N2O generation
nylon46bioreductHierarchist <- nylonComparisons(results,comparisonFile=comparisonFiles[["ComparisonsHierarchist"]],
                                         parameters,MolPar=allParameters,linkages,flowData,
                                         backgroundData,CFs,data.dir,mainNylon="Nylon46Reduct",refNylon = "Nylon66Reduct",perspective="Hierarchist",
                                         main="Putrescine",reactant="AdipicAcid_average",reference=c("AdipicAcid_average","HMDA_Glo"))
nylon46bioreductNormHierarchist <- normalisation(data.dir="Data/ReCiPe",nylon46bioreductHierarchist@aggregated[nylon46bioreductHierarchist@aggregated$mol=="Collated",],type="nylon")

#Normalise each of the nylon scenarios
allRes <- list()
allResNorm <- list()
for(ii in 1:length(perspectives)) {
  nylon46Restemp <- nylonComparisons(results,comparisonFile=comparisonFiles[[paste0("Comparisons",perspectives[ii])]],parameters,MolPar=allParameters,linkages,flowData,
                                 backgroundData,CFs,data.dir,mainNylon="Nylon46",refNylon="Nylon66",perspective=perspectives[ii],
                                 main="Putrescine",reactant="AdipicAcid_RoW",reference=c("AdipicAcid_RoW","HMDA_Glo"))
  nylon46ResNorm <- normalisation(data.dir="Data/ReCiPe",nylon46Restemp@aggregated[nylon46Restemp@aggregated$mol=="Collated",],type="nylon")
  assign(paste0("nylon46Res",perspectives[ii]),nylon46Restemp)
  assign(paste0("nylon46ResNorm",perspectives[ii]),nylon46ResNorm)
  
  nylon410Restemp <- nylonComparisons(results,comparisonFile=comparisonFiles[[paste0("Comparisons",perspectives[ii])]],
                                      parameters,MolPar=allParameters,linkages,flowData,
                                  backgroundData,CFs,data.dir,mainNylon="Nylon410",refNylon="Nylon66",perspective=perspectives[ii],
                                  main="Putrescine",reactant="SebacicAcid_average",reference=c("AdipicAcid_RoW","HMDA_Glo"))
  nylon410ResNorm <- normalisation(data.dir="Data/ReCiPe",nylon410Restemp@aggregated[nylon410Restemp@aggregated$mol=="Collated",],type="nylon")
  assign(paste0("nylon410Res",perspectives[ii]),nylon410Restemp)
  assign(paste0("nylon410ResNorm",perspectives[ii]),nylon410ResNorm)
  
  nylon510Restemp <- nylonComparisons(results,comparisonFile=comparisonFiles[[paste0("Comparisons",perspectives[ii])]],
                                  parameters,MolPar=allParameters,linkages,flowData,
                                  backgroundData,CFs,data.dir,mainNylon="Nylon510",refNylon="Nylon66",perspective=perspectives[ii],
                                  main="Cadaverine",reactant="SebacicAcid_average",reference=c("AdipicAcid_RoW","HMDA_Glo"))
  nylon510ResNorm <- normalisation(data.dir="Data/ReCiPe",nylon510Restemp@aggregated[nylon510Restemp@aggregated$mol=="Collated",],type="nylon")
  assign(paste0("nylon510Res",perspectives[ii]),nylon510Restemp)
  assign(paste0("nylon510ResNorm",perspectives[ii]),nylon510ResNorm)
  
  allRes[[ii]] <- rbind(nylon46Restemp@aggregated,nylon410Restemp@aggregated[nylon410Restemp@aggregated$key == "main",],
                        nylon510Restemp@aggregated[nylon510Restemp@aggregated$key == "main",])
  
  allResNorm[[ii]] <- rbind(nylon46ResNorm,nylon410ResNorm[nylon410ResNorm$key == "main",],nylon510ResNorm[nylon510ResNorm$key == "main",])
  
}
allNylon <- allRes %>% Reduce(function(df1,df2) rbind(df1,df2),.)
allNylonNorm <- allResNorm %>% Reduce(function(df1,df2) rbind(df1,df2),.)
allNylonNorm$Scenario <- paste0(allNylonNorm$nylon,allNylonNorm$Perspective)

#Variability test
variabilityCadaverine <- sensTest(data.dir,baseResult=results,sensTestPars,linkages,parameters,CFs,target="Cadaverine",perspective="Hierarchist",type="MinMax",MCruns=MCruns,MultiStart=MultiStart)
write.csv(variabilityCadaverine,file=file.path(output.dir,paste0("VariabilityCadaverine.csv")))
variabilityPutrescine <- sensTest(data.dir,baseResult=results,sensTestPars,linkages,parameters,CFs,target="Putrescine",perspective="Hierarchist",type="MinMax",MCruns=MCruns,MultiStart=MultiStart)
write.csv(variabilityPutrescine,file=file.path(output.dir,paste0("VariabilityPutrescine.csv")))


save(results,resultsNormalised,variabilityPutrescine,variabilityCadaverine,
     nylon46ResHierarchist,nylon410ResHierarchist,nylon510ResHierarchist,
     nylon46ResEgalitarian,nylon410ResEgalitarian,nylon510ResEgalitarian,
     nylon46ResIndividualist,nylon410ResIndividualist,nylon510ResIndividualist,
     nylon46bioreductHierarchist,nylon46bioreductNormHierarchist,
     allNylon,allNylonNorm,parameters,
     file=file.path(output.dir,paste0(overallName,".RData")))
