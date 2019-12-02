######Set-up#####
library(ggplot2)
library(xlsx)
library(stringr)
library(magrittr)
library(tidyr)

#Source my own themes
source("Scripts/PublicationPlots.R")
source("Scripts/SustAssessRSourceCode.R")

#set data directory
input.dir <- "Data/Results"
input.file <- "CadaverinePutrescine.Rdata"
#Set and create output directory
output.dir <- "ModelOutputs"
dir.create(file.path(output.dir),showWarnings = FALSE)

#Load data
load(file.path(input.dir,input.file))

#####Summarise all flows#####
targets=c("Putrescine","Cadaverine");perspective=c("Hierarchist","Egalitarian","Individualist")

resultNames <- names(results@RawLCAresult[grepl(paste(targets,collapse="|"),names(results@RawLCAresult)) &
                                          grepl(paste(perspective,collapse="|"),names(results@RawLCAresult))])
impactCatUnits <- c("kg CO2 eq.","kg PM2.5 eq.","kg oil eq.","m3","kg 1,4 DB eq.","kg P eq.","kg 1,4-DB eq.",
                    "kg 1,4-DB eq.","Bq C-60 eq. to air","Annual crop eq.?y","kg 1,4-DB eq.","kg N eq.","kg Cu eq.",
                    "kg NOx eq.","kg NOx eq.","kg CFC-11 eq.","kg SO2 eq.","kg 1,4-DB eq.")
flowresults <- list()
for(ii in 1:length(resultNames)) {
  #Flow results
  tempRes <- results@RawLCAresult[[resultNames[ii]]]
  inputNames <- names(tempRes@inputs)
  inputList <- list()
  outputNames <- names(tempRes@outputs)
  outputList <- list()
  #Set up indices for cleaning
  ModuleNameIndex <- c(feedstockUSLigno = "Ligno Sugar Production",distillation2="Distillation",distillation1="Stripper",
                   extractor="Extractor",phAdjustment="pH Adjustment",centrifuge="Centrifuge",fermenter="Fermenter",
                   dryer="Dryer",evaporatorBiomass="Evaporator (Biomass)",evaporatorCake="Evaporator (Cake)",combustor="Combustor")
  for(jj in 1:length(tempRes@inputs)) {
    tempStats <- extractStats(tempRes@inputs[[jj]])
    #Put in the necessary info
    tempStats <- cbind(rep(resultNames[ii],nrow(tempStats)),rep(inputNames[jj],nrow(tempStats)),
                       rep("Input",nrow(tempStats)),row.names(tempStats),rep("kg",nrow(tempStats)),
                       tempStats[2:ncol(tempStats)],stringsAsFactors=FALSE)
    colnames(tempStats) <- c("Result Name","Module","Flow Type","Flow","Units","Mean","Min","Max","CI95low","CI95up",
                  "Variance","Standard Deviation","Variability","Skewness","Kurtosis")
    inputList[[jj]] <- tempStats
  }
  inputs <- inputList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  for(jj in 1:length(tempRes@outputs)) {
    tempStats <- extractStats(tempRes@outputs[[jj]])
    #Put in the necessary info
    tempStats <- cbind(rep(resultNames[ii],nrow(tempStats)),rep(inputNames[jj],nrow(tempStats)),
                       rep("Output",nrow(tempStats)),row.names(tempStats),rep("kg",nrow(tempStats)),
                       tempStats[2:ncol(tempStats)],stringsAsFactors=FALSE)
    colnames(tempStats) <- c("Result Name","Module","Flow Type","Flow","Units","Mean","Min","Max","CI95low","CI95up",
                             "Variance","Standard Deviation","Variability","Skewness","Kurtosis")
    outputList[[jj]] <- tempStats
  }  
  outputs <- outputList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  inputoutput <- rbind(inputs,outputs)
  #Correct units
  inputoutput$Units[inputoutput$Flow == "Electricity"] <- "kWh"
  inputoutput$Units[inputoutput$Flow == "Steam"] <- "MJ"

  #Remove product if not for distillation2
  inputoutput <- inputoutput[!(inputoutput$Module != "distillation2" & inputoutput$Flow == "Product"),]
  
  #Replace the module names
  inputoutput$Module <- ModuleNameIndex[inputoutput$Module]
  
  
  #Write to file
  if(ii==1){
    write.xlsx(inputoutput, file=file.path(output.dir,paste0("ProcessFlowSummaries",".xlsx")),
               sheetName=resultNames[ii], append=FALSE, row.names = FALSE)
  } else {
    write.xlsx(inputoutput, file=file.path(output.dir,paste0("ProcessFlowSummaries",".xlsx")),
               sheetName=resultNames[ii], append=TRUE, row.names = FALSE)
  }
  
  #Aggregated LCI
  testListInputs <- list(RawMaterials = c("Sugar","Salt","Nitrogen","Water","Butanol","NaOH"),
                         Utilities = c("Steam","Electricity"),
                         LignoSystemInput = c("Lime","CornStover","SulfuricAcid","CornSteepLiqour","NaOH","Ammonia","Glucose","AmmoniumSulfate","SoybeanOil","SO2"),
                         other = c("CombustorUsage","FermPlant","EmbodiedCarbon"))
                         
                   
  testListOutputs <- list(Waste=c("Waste"),
                          Emissions = c("CO2Emit","COEmit","SO2Emit","NO2Emit","CH4Emit"),
                          UtilityCredit="Electricity")
  allResList <- list()
  #Inputs
  for(kk in 1:length(testListInputs)) {
    tempTest <- testListInputs[[kk]]
    tempResList <- list()
    for(jj in 1:length(tempTest)){
      totalVol <- 0
      #Add in usage for every module
      for(mm in which(sapply(tempRes@inputs,function(x) tempTest[jj] %in% names(x)))) {
        totalVol <- totalVol + tempRes@inputs[[mm]][[tempTest[jj]]]
      }
      #Minus if it's used-up elsewhere
      for(mm in which(sapply(tempRes@outputs,function(x) tempTest[jj] %in% names(x)))) {
        totalVol <- totalVol - tempRes@outputs[[mm]][[tempTest[jj]]]
      }
      #Do corrections for steam
      if(tempTest[jj] == "Steam") {
        #Do not allow steam to be in credit
        totalVol[totalVol < 0] <- 0
      }
      #Replace very small numbers with 0
      totalVol[totalVol < 1e-10] <- 0
      tempResList[[tempTest[jj]]] <- totalVol
    }
    allResList[[names(testListInputs)[kk]]] <- cbind(data.frame(Scenario=resultNames[ii],InputOutput="Input",Type=names(testListInputs)[kk],Flow=testListInputs[[kk]],Units="kg",stringsAsFactors = FALSE),extractStats(tempResList))
  }
  
  #Ouputs
  for(kk in 1:length(testListOutputs)) {
    tempTest <- testListOutputs[[kk]]
    tempResList <- list()
    for(jj in 1:length(tempTest)){
      totalVol <- 0
      #Add in usage for every module
      for(mm in which(sapply(tempRes@outputs,function(x) tempTest[jj] %in% names(x)))) {
        totalVol <- totalVol + tempRes@outputs[[mm]][[tempTest[jj]]]
      }
      tempResList[[tempTest[jj]]] <- totalVol
      #Replace very small numbers with 0
      totalVol[totalVol < 1e-10] <- 0
    }
    allResList[[names(testListOutputs)[kk]]] <- cbind(data.frame(Scenario=resultNames[ii],InputOutput="Output",Type=names(testListOutputs)[kk],Flow=testListOutputs[[kk]],Units="kg",stringsAsFactors = FALSE),extractStats(tempResList))
  }
  LCI <- allResList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  LCI$indicator <- NULL
  LCI$Units[LCI$Flow == "Electricity"] <- "kWh"
  LCI$Units[LCI$Flow == "Steam"] <- "MJ"
  
  if(ii==1){
    write.xlsx(LCI, file=file.path(output.dir,paste0("AggregatedLCI",".xlsx")),
               sheetName=resultNames[ii], append=FALSE, row.names = FALSE)
  } else {
    write.xlsx(LCI, file=file.path(output.dir,paste0("AggregatedLCI",".xlsx")),
               sheetName=resultNames[ii], append=TRUE, row.names = FALSE)
  }
  
  #Hotspot results
  hotspotResults <- results@LCAHotspot[[resultNames[ii]]]
  hotspotResults$Indicator <- row.names(hotspotResults)
  hotspotResults <- hotspotResults[hotspotResults$Indicator != "Climate change, incl biogenic carbon",]
  
  hotspotResults$Units <- impactCatUnits
  #Reorder hotspot results
  hotspotResults <- hotspotResults[,c((ncol(hotspotResults)-1):ncol(hotspotResults),1:(ncol(hotspotResults)-2))]
  hotspotResults <- hotspotResults[,c(TRUE,TRUE,abs(colSums(hotspotResults[3:ncol(hotspotResults)])) > 0)]
  hotspotResults$Indicator[hotspotResults$Indicator == "Climate change, default, excl biogenic carbon"] <- "Climate change (with credit for embodied carbon)"
  
  if(ii==1){
    write.xlsx(hotspotResults, file=file.path(output.dir,paste0("HotspotSummaries",".xlsx")),
               sheetName=resultNames[ii], append=FALSE, row.names = FALSE)
  } else {
    write.xlsx(hotspotResults, file=file.path(output.dir,paste0("HotspotSummaries",".xlsx")),
               sheetName=resultNames[ii], append=TRUE, row.names = FALSE)
  }
}

#Export all LCIA results
#For the monomers
col.order <- c("Monomer","Perspective","indicator","Units","mean","min","max","CI95low","CI95up",
               "var","sd","variability",
               "skewness","kurtosis")
col.names.manual <- c("Monomer","Perspective","Indicator","Units","Mean","Min","Max","CI95low","CI95up",
                      "Variance","Standard deviation","Variability",
                      "Skewness","Kurtosis")

allMonExport <- results@LCAresultsDF[results@LCAresultsDF$indicator != "Climate change, incl biogenic carbon",]
allMonExport$indicator[allMonExport$indicator == "Climate change, default, excl biogenic carbon"] <- "Climate change (with credit for embodied carbon)"
allMonExport$Units <- rep(impactCatUnits,nrow(allMonExport)/length(impactCatUnits))
allMonExport$Monomer <- str_remove(allMonExport$Scenario,"Hierarchist|Individualist|Egalitarian")
allMonExport$Perspective <- str_remove(allMonExport$Scenario,"Cadaverine|Putrescine")
allMonExport <- allMonExport[,col.order]
colnames(allMonExport) <- col.names.manual
write.csv(allMonExport,file=file.path(output.dir,paste0("MonomerLCIAResults",".csv")),row.names=FALSE)
write.xlsx(allMonExport,file=file.path(output.dir,paste0("MonomerLCIAResults",".xlsx")),append=FALSE, row.names = FALSE,sheetName = "MonomerLCIAResults")

#Export nylon results
col.order <- c("Scenario","Perspective","nylon","mol","indicator","Units","mean","min","max","CI95low","CI95up",
               "var","sd","variability",
               "skewness","kurtosis")
col.names.manual <- c("Scenario","Perspective","Nylon","Component","Indicator","Units","Mean","Min","Max","CI95low","CI95up",
                      "Variance","Standard deviation","Variability",
                      "Skewness","Kurtosis")
allNylonExtra <- rbind(allNylon[allNylon$Perspective=="Hierarchist",],nylon46bioreductHierarchist@aggregated)
allNylonExtra <- allNylonExtra[allNylonExtra$nylon %in% c("Nylon410","Nylon510", "Nylon46Reduct","Nylon66Reduct"),]
allNylonExport <- allNylonExtra[allNylonExtra$indicator != "Climate change, incl biogenic carbon",]
allNylonExport$Units <- rep(impactCatUnits,nrow(allNylonExport)/length(impactCatUnits))
allNylonExport <- allNylonExport[,col.order]
allNylonExport$indicator[allNylonExport$indicator == "Climate change, default, excl biogenic carbon"] <- "Climate change (with credit for embodied carbon)"

colnames(allNylonExport) <- col.names.manual
allNylonExport <- mutate_if(allNylonExport,is.numeric, ~replace(., is.na(.), 0))
write.csv(allNylonExport,file=file.path(output.dir,paste0("AllNylonLCIAResults",".csv")),row.names = FALSE)
write.xlsx(allNylonExport,file=file.path(output.dir,paste0("AllNylonLCIAResults",".xlsx")),row.names = FALSE,append=FALSE,sheetName = "AllNylonLCIAResults")

#Export cost results
costBreakdown <- results@CostBreakdownDF
costBreakdown$Units <- "$"
costBreakdown$indicator <- costBreakdown$Area
costBreakdown$Area <- NULL

col.order <- c("Scenario","indicator","Units","mean","min","max","CI95low","CI95up",
               "var","sd","variability",
               "skewness","kurtosis")
col.names.manual <- c("Scenario","Cost Area","Units","Mean","Min","Max","CI95low","CI95up",
                      "Variance","Standard deviation","Variability",
                      "Skewness","Kurtosis")
MSPdf <- results@MSP_DF
MSPdf$Units <- "$/kg"
costResults <- rbind(costBreakdown,MSPdf)
costResults$Scenario <- str_remove(costResults$Scenario,"Hierarchist")
costResults <- costResults[,col.order]
colnames(costResults) <- col.names.manual

write.csv(costResults,file=file.path(output.dir,paste0("CostBreakdown",".csv")),row.names = FALSE)
write.xlsx(costResults,file=file.path(output.dir,paste0("CostBreakdown",".xlsx")),row.names = FALSE,append=FALSE,sheetName = "CostBreakdown")

#Export Sensitivity results
col.order <- c("Scenario","indicator","par","type","mean","min","max","CI95low","CI95up",
               "var","sd","variability",
               "skewness","kurtosis")
col.names.manual <- c("Scenario","Indicator","Parameter","Parameter Type","Mean","Min","Max","CI95low","CI95up",
               "Variance","Standard deviation","Variability",
               "Skewness","Kurtosis")
variabilityPutrescine4Export <- variabilityPutrescine[variabilityPutrescine$mean != 0 & variabilityPutrescine$indicator != "Climate change, incl biogenic carbon",
                                col.order]
variabilityPutrescine4Export$indicator[variabilityPutrescine4Export$indicator == "Climate change, default, excl biogenic carbon"] <- "Climate change (with credit for embodied carbon)"
colnames(variabilityPutrescine4Export) <- col.names.manual
write.xlsx(variabilityPutrescine4Export, file=file.path(output.dir,paste0("Sensitivity",".xlsx")),
           sheetName="Putrescine", append=FALSE, row.names = FALSE)


variabilityCadaverine4Export <- variabilityCadaverine[variabilityCadaverine$mean != 0 & variabilityCadaverine$indicator != "Climate change, incl biogenic carbon",
                                                      col.order]
variabilityCadaverine4Export$indicator[variabilityCadaverine4Export$indicator == "Climate change, default, excl biogenic carbon"] <- "Climate change (with credit for embodied carbon)"
colnames(variabilityCadaverine4Export) <- col.names.manual
write.xlsx(variabilityCadaverine4Export, file=file.path(output.dir,paste0("Sensitivity",".xlsx")),
           sheetName="Cadaverine", append=TRUE, row.names = FALSE)

#####Nylon usage comparisons#####
#Just for climate change, hierarchist
dataset <- rbind(allNylon[!(allNylon$nylon %in% c("Nylon66","Nylon46")) & allNylon$Perspective == "Hierarchist" & allNylon$indicator == "Climate change, default, excl biogenic carbon",],
                 nylon46bioreductHierarchist@aggregated[nylon46bioreductHierarchist@aggregated$indicator == "Climate change, default, excl biogenic carbon",])
dataset <- dataset[dataset$mol != "Collated",]

nylonNames <- c("Nylon 66\n(Fossil)","Nylon 46\n(Part Bio)","Nylon 410\n(Full Bio)","Nylon 510\n(Full Bio)")
dataset$nylon <- factor(dataset$nylon,levels=c("Nylon66Reduct","Nylon46Reduct","Nylon410","Nylon510"))
dataset$mol <- factor(dataset$mol,levels=c("Polymerisation","AdipicAcid_average","AdipicAcid_Bio","HMDA_Glo","SebacicAcid_average","Putrescine","Cadaverine"))
legendLabs=c("Polymerisation","Adipic Acid (Fossil)","HMDA (Fossil)","SebacicAcid (Bio)","Putrescine (Bio)","Cadaverine (Bio)")

CIRes <- data.frame(nylon=rep(levels(dataset$nylon),each=1),
                    mean=dataset$mean[dataset$mol=="Polymerisation"],
                    CI95low=unlist(lapply(levels(dataset$nylon),function(x) sum(dataset$CI95low[dataset$nylon==x]))),
                    CI95up=unlist(lapply(levels(dataset$nylon),function(x) sum(dataset$CI95up[dataset$nylon==x]))))
CIRes$mol <- "Polymerisation"


ggplot(data = dataset,aes(x = nylon, y = mean, fill = mol)) + 
  geom_bar(stat="identity") +
  geom_errorbar(data=CIRes,aes(ymin=CI95low,ymax=CI95up),color="black",position=position_dodge(0.9),width=.4,size=2) +
  theme_Publication() +
  scale_fill_Publication(labels=legendLabs) +
  labs(x="\nNylon Scenario", y="Climate Change\n(kg CO2 eq./kg Nylon)") +
  scale_x_discrete(labels=nylonNames) +
  theme(axis.title = element_text(size=25),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=18),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=25,face="bold"),
        legend.text = element_text(size=20),
        panel.spacing = unit(1, "lines"),strip.background = element_blank(),strip.text = element_text(size=25,face="bold")) +
  guides(fill=guide_legend(
    title="Component",
    keywidth=0.3,
    keyheight=0.3,
    default.unit="inch"))
ggsave(file.path(output.dir,"NylonUsageClimateChangeHierarchist.pdf"),width = 14, height = 7.5,dpi=500)
ggsave(file.path(output.dir,"NylonUsagesClimateChangeHierarchist.png"),width = 14, height = 7.5,dpi=500)



#####Plot comparative results - different nylons#####
#Declare function for normalisation
maxNorm <- function(allData) {
  maxVal <- c()
  indicators <- unique(allData$indicator)
  for(ii in 1:length(indicators)) {
    maxVal[ii] <- max(allData$mean[allData$indicator == indicators[ii]])
  }
  normFactors <- rep(maxVal,length(unique(allData$Scenario)))
  allDataNorm <- allData
  allDataNorm$mean <- allDataNorm$mean/normFactors*100
  allDataNorm$min <- allDataNorm$min/normFactors*100
  allDataNorm$max <- allDataNorm$max/normFactors*100
  allDataNorm$CI95up <- allDataNorm$CI95up/normFactors*100
  allDataNorm$CI95low <- allDataNorm$CI95low/normFactors*100
  return(allDataNorm)
}

impactCat <- c("Climate change, default, excl biogenic carbon", 
               "Fossil depletion","Metal depletion","Freshwater Consumption","Freshwater ecotoxicity","Freshwater Eutrophication",
               "Human toxicity, cancer","Human toxicity, non-cancer","Fine Particulate Matter Formation","Ionizing Radiation","Land use","Marine ecotoxicity",
               "Marine Eutrophication", 
               "Photochemical Ozone Formation, Ecosystems","Photochemical Ozone Formation, Human Health","Stratospheric Ozone Depletion",
               "Terrestrial Acidification","Terrestrial ecotoxicity")

impactCatNames <- c("CC","FD","MD","FC","FE","FET","HT\n(C)","HT\n(NC)","PM","IR","LU","ME","MET","POF\n(Eco)",
                    "POF\n(HH)","SOD","TA","TE")

perspectives <- unique(allNylon$Perspective)
for(ii in 1:length(perspectives)){
  dataset <- allNylon[allNylon$mol == "Collated"&allNylon$Perspective==perspectives[ii],]
  datasetNorm <- maxNorm(dataset)
  nylonScens <- c("Nylon66","Nylon46","Nylon410","Nylon510")
  datasetNorm$nylon <- factor(datasetNorm$nylon,levels=nylonScens)
  datasetNormSel <- datasetNorm[datasetNorm$indicator %in% impactCat,]
  ggplot(data = datasetNormSel,aes(x = indicator, y = mean, fill = nylon)) + 
    geom_bar(stat="identity", position = "dodge") +
    geom_errorbar(aes(ymin=CI95low,ymax=CI95up),color="black",position=position_dodge(.9),width=.4,size=2) +
    theme_Publication() +
    scale_fill_Publication(name="Nylons",labels=nylonScens) +
    labs(x="Impact Category (ReCiPe)", y="Result\n(% maximum for\nimpact category)") +
    ggtitle(paste("Comparative Impacts",unique(dataset$geog),unique(dataset$nylon),"\n")) +
    scale_x_discrete(labels=impactCatNames) +
    theme(axis.title = element_text(size=28),
          axis.text.y = element_text(size=35),
          axis.text.x = element_text(size=25),
          plot.title = element_text(hjust = 0.5, size = 45),
          legend.title = element_text(size=28,face="bold"),
          legend.text = element_text(size=28)) +
    guides(fill=guide_legend(
      keywidth=0.3,
      keyheight=0.3,
      default.unit="inch"))
  ggsave(file.path(output.dir,paste0("AllImpactCategories",perspectives[ii],".pdf")),width = 18, height = 8.5,dpi=500)
  ggsave(file.path(output.dir,paste0("AllImpactCategories",perspectives[ii],".png")),width = 18, height = 8.5,dpi=500)

}

#Plot comparative results, just hierarchist, with adipic acid reduced impact
allNylonExtra <- rbind(allNylon[allNylon$Perspective=="Hierarchist"& allNylon$mol=="Collated",],nylon46bioreductHierarchist@aggregated[nylon46bioreductHierarchist@aggregated$mol=="Collated",])
dataset <- allNylonExtra[allNylonExtra$nylon %in% c("Nylon410","Nylon510", "Nylon46Reduct","Nylon66Reduct"),]
datasetNorm <- maxNorm(dataset)
nylonScens <- c("Nylon66Reduct","Nylon46Reduct","Nylon410","Nylon510")
nylonScenNames <- c("Nylon 66","Nylon 46","Nylon 410","Nylon 510")
datasetNorm$nylon <- factor(datasetNorm$nylon,levels=nylonScens)
datasetNormSel <- datasetNorm[datasetNorm$indicator %in% impactCat,]
ggplot(data = datasetNormSel,aes(x = indicator, y = mean, fill = nylon)) + 
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin=CI95low,ymax=CI95up),color="black",position=position_dodge(.9),width=.4,size=2) +
  theme_Publication() +
  scale_fill_Publication(name="Nylons",labels=nylonScenNames) +
  labs(x="Impact Category (ReCiPe)", y="Result (% max. for\nimpact category)") +
  ggtitle(paste("Comparative Impacts","\n")) +
  scale_x_discrete(labels=impactCatNames) +
  theme(axis.title = element_text(size=28),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=25),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=28,face="bold"),
        legend.text = element_text(size=28)) +
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.3,
    default.unit="inch"))
ggsave(file.path(output.dir,paste0("AllImpactCategoriesHierarchistReduct",".pdf")),width = 16, height = 7.5,dpi=500)
ggsave(file.path(output.dir,paste0("AllImpactCategoriesHierarchistReduct",".png")),width = 16, height = 7.5,dpi=500)

#####Plot hotspots#####
targets="Putrescine";perspective="Hierarchist"

impactCat <- c("Climate change, default, excl biogenic carbon", 
               "Fossil depletion","Metal depletion","Freshwater Consumption","Freshwater ecotoxicity","Freshwater Eutrophication",
               "Human toxicity, cancer","Human toxicity, non-cancer","Fine Particulate Matter Formation","Ionizing Radiation","Land use","Marine ecotoxicity",
               "Marine Eutrophication", 
               "Photochemical Ozone Formation, Ecosystems","Photochemical Ozone Formation, Human Health","Stratospheric Ozone Depletion",
               "Terrestrial Acidification","Terrestrial ecotoxicity")

impactCatNames <- c("CC","FD","MD","FC","FE","FET","HT\n(C)","HT\n(NC)","PM","IR","LU","ME","MET","POF\n(Eco)",
                    "POF\n(HH)","SOD","TA","TE")

#Select the necessary results to use
resultsUsed <- results@LCAHotspot[which(grepl(paste(targets,collapse="|"),names(results@LCAHotspot)) & 
                                          grepl(paste(perspective,collapse="|"),names(results@LCAHotspot)))]

#Extract and normalise results
allDataNorm <- resultsUsed[[1]]
allDataNorm <- allDataNorm/rowSums(allDataNorm)*100
#Remove any with nothing
allDataNorm <- allDataNorm[,which(colSums(allDataNorm)!=0)]
allDataNorm$Indicator <- rownames(allDataNorm)
#And gather it up
allDataGathered <- gather(allDataNorm,key="UnitProcess",value="Value",1:(ncol(allDataNorm)-1))
tempList <- list()
indicators <- unique(allDataGathered$Indicator)
#Remove and aggregate any that are less than 5%
for(ii in 1:length(indicators)) {
  temp <- allDataGathered[allDataGathered$Indicator == indicators[ii],]
  tempFiltered <- temp[abs(temp$Value) >= 5,]
  tempFiltered[nrow(tempFiltered)+1,] <- c(indicators[ii],"Others",sum(temp$Value[abs(temp$Value) < 5]))
  tempList[[ii]] <- tempFiltered
}
resultsDF <- tempList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
resultsDF$Value <- as.numeric(resultsDF$Value)
resultsDF$UnitProcess <- factor(resultsDF$UnitProcess,levels=c(unique(resultsDF$UnitProcess)[unique(resultsDF$UnitProcess)!="Others"],"Others"))
resultsDFPrep <- resultsDF[resultsDF$Indicator %in% impactCat,]
resultsDFPrep$Indicator <- factor(resultsDFPrep$Indicator,levels=impactCat)
ggplot(data = resultsDFPrep,aes(x = Indicator, y = Value, fill = UnitProcess)) + 
  geom_bar(stat="identity") +
  theme_Publication() +
  scale_fill_Publication() +
  labs(x="Impact Category (ReCiPe)", y="Result\n(% contribution)") +
  ggtitle(paste("Hotspots","\n")) +
  scale_x_discrete(labels=impactCatNames) +
  theme(axis.title = element_text(size=28),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=24),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=28,face="bold"),
        legend.text = element_text(size=24)) +
  guides(fill=guide_legend(
    keywidth=0.3,
    keyheight=0.3,
    default.unit="inch"))
ggsave(file.path(output.dir,paste0("Putrescine","Hotspots",".pdf")),width = 18, height = 8.5,dpi=500)
ggsave(file.path(output.dir,paste0("Putrescine","Hotspots",".png")),width = 18, height = 8.5,dpi=500)

#####Variability results#####
#First for flow parameters
impactCat <- c("Climate change, default, excl biogenic carbon", 
               "Fossil depletion","Metal depletion","Freshwater Consumption","Freshwater ecotoxicity","Freshwater Eutrophication",
               "Human toxicity, cancer","Human toxicity, non-cancer","Fine Particulate Matter Formation","Ionizing Radiation","Land use","Marine ecotoxicity",
               "Marine Eutrophication", 
               "Photochemical Ozone Formation, Ecosystems","Photochemical Ozone Formation, Human Health","Stratospheric Ozone Depletion",
               "Terrestrial Acidification","Terrestrial ecotoxicity")

impactCatNames <- c("CC","FD","MD","FC","FE","FET","HT\n(C)","HT\n(NC)","PM","IR","LU","ME","MET","POF\n(Eco)",
                    "POF\n(HH)","SOD","TA","TE")

allData <- variabilityPutrescine
allData <- allData[allData$type == "flow",]
tempList <- list()
indicators <- unique(allData$indicator)
#Remove and aggregate any that are less than 5%
for(ii in 1:length(indicators)) {
  temp <- allData[allData$indicator == indicators[ii],]
  tempFiltered <- temp[abs(temp$mean) >= 0.05,]
  tempFiltered[nrow(tempFiltered)+1,] <- c(indicators[ii],
                                           sum(temp$mean[abs(temp$mean) < 0.05]),sum(temp$min[abs(temp$mean) < 0.05]),
                                           sum(temp$max[abs(temp$mean) < 0.05]),sum(temp$var[abs(temp$mean)< 0.05]),
                                           sum(temp$variability[abs(temp$mean) < 0.05]),sum(temp$sd[abs(temp$mean) < 0.05]),
                                           sum(temp$skewness[abs(temp$mean) < 0.05]),sum(temp$kurtosis[abs(temp$mean) < 0.05]),
                                           sum(temp$CI95up[abs(temp$mean) < 0.05]),sum(temp$CI95low[abs(temp$mean) < 0.05]),
                                           temp$Scenario[1],"Others",unique(allData$type))
  tempList[[ii]] <- tempFiltered
}
resultsDF <- tempList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
resultsDF$mean <- as.numeric(resultsDF$mean)
resultsDFPrep <- resultsDF[resultsDF$indicator %in% impactCat,]
resultsDFPrep$indicator <- factor(resultsDFPrep$indicator,levels=impactCat)
parsUsed <- c("MaxYield","OrgSpeed","Switch","DSPs","energySource","NitrogenSource","geogs","WasteScenario","Others")
resultsDFPrep$par <- factor(resultsDFPrep$par,levels=parsUsed)
names <- c("Yield on\nSugar","Organism\nSpeed","Switch \ntime",
                    "Downstream Processing\nEfficiency", "Energy Source\n(Grid or biomass)",
                    "Nitrogen\nSource","Location &\nFeedstock",
                    "Integrated vs\nNon-Integrated",
                    "Other\nParameters")
resultsDFPrep$CI95lowNew <- as.numeric(resultsDFPrep$CI95low)
resultsDFPrep$CI95upNew <- as.numeric(resultsDFPrep$CI95up)
resultsDFPrep$test <- rep(0.1,nrow(resultsDFPrep))

resultsDFPrep$mean <- resultsDFPrep$mean*100
ggplot(data=resultsDFPrep,aes(x=indicator,y=mean,fill=par)) +
  geom_bar(stat="identity",position="dodge") +
  theme_Publication() +
  scale_fill_Publication(name="Parameter",labels=names) +
  labs(x="Impact Category (ReCiPe)", y="Sensitivity (%)") +
  ggtitle(paste("Sensitivity Analysis\n")) +
  scale_x_discrete(labels = impactCatNames) +
  theme(axis.title = element_text(size=28),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=24),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=26,face="bold"),
        legend.text = element_text(size=20))+
  guides(fill=guide_legend(
    keywidth=0.6,
    keyheight=0.6,
    default.unit="inch"))
ggsave(file.path(output.dir,paste0("FlowSensitivityPutrescine","Comparisons",".pdf")),width = 18, height = 8.5,dpi=500)
ggsave(file.path(output.dir,paste0("FlowSensitivityPutrescine","Comparisons",".png")),width = 18, height = 8.5,dpi=500)

#Second for process parameters
allData <- variabilityPutrescine
allData <- allData[allData$type == "process",]
tempList <- list()
indicators <- unique(allData$indicator)
#Remove and aggregate any that are less than 5%
for(ii in 1:length(indicators)) {
  temp <- allData[allData$indicator == indicators[ii],]
  tempFiltered <- temp[abs(temp$mean) >= 0.05,]
  tempFiltered[nrow(tempFiltered)+1,] <- c(indicators[ii],
                                           sum(temp$mean[abs(temp$mean) < 0.05]),sum(temp$min[abs(temp$mean) < 0.05]),
                                           sum(temp$max[abs(temp$mean) < 0.05]),sum(temp$var[abs(temp$mean) < 0.05]),
                                           sum(temp$variability[abs(temp$mean) < 0.05]),sum(temp$sd[abs(temp$mean) < 0.05]),
                                           sum(temp$skewness[abs(temp$mean) < 0.05]),sum(temp$kurtosis[abs(temp$mean) < 0.05]),
                                           sum(temp$CI95up[abs(temp$mean) < 0.05]),sum(temp$CI95low[abs(temp$mean) < 0.05]),
                                           temp$Scenario[1],"Others",unique(allData$type))
  tempList[[ii]] <- tempFiltered
}
resultsDF <- tempList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
resultsDF$mean <- as.numeric(resultsDF$mean)
resultsDFPrep <- resultsDF[resultsDF$indicator %in% impactCat,]
resultsDFPrep$indicator <- factor(resultsDFPrep$indicator,levels=impactCat)

names <- str_wrap(c("Biorefinery Size","Feedstock Input","NaOH Input",
                    "Solvent Input", "Waste Generation","Water Input","Water Evaporated",
                    "Other Parameters"),20)

resultsDFPrep$mean <- resultsDFPrep$mean*100

ggplot(data=resultsDFPrep,aes(x=indicator,y=mean,fill=par)) +
  geom_bar(stat="identity",position="dodge") +
  theme_Publication() +
  scale_fill_Publication(name="Parameter") +
  labs(x="Impact Category (ReCiPe)", y="\nSensitivity (%)") +
  ggtitle(paste("Sensitivity Analysis")) +
  scale_x_discrete(labels = impactCatNames) +
  theme(axis.title = element_text(size=25),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=24),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=25,face="bold"),
        legend.text = element_text(size=20))+
  guides(fill=guide_legend(
    keywidth=0.5,
    keyheight=0.5,
    default.unit="inch"))
ggsave(file.path(output.dir,paste0("ProcessSensitivityPutrescine","Comparisons",".pdf")),width = 18, height = 9,dpi=500)
ggsave(file.path(output.dir,paste0("ProcessSensitivityPutrescine","Comparisons",".png")),width = 18, height = 9,dpi=500)

#Cost variability
pars <- c("MaxYield","OrgSpeed","Switch","DSPs","geogs","BioethanolPlant","Sugar","Nitrogen","discountRate","times")

parNames <- c("Yield on\nSugar","Organism\nSpeed","Switch\ntime",
           "Downstream\nProcessing\nEfficiency","Location\nand\nFeedstock",
           "Biorefinery\nBase Cost","Sugar\nPrice","Ligno Sugar\nPrice","Nitrogen\nPrice","Discount\nRate","Time\nspan")

MSPvariability <- variabilityPutrescine[variabilityPutrescine$indicator=="MSP" & variabilityPutrescine$par != "overall",]
MSPvariability <- MSPvariability[MSPvariability$mean > 0.05,]
MSPvariability$par <- factor(MSPvariability$par,levels = pars)
MSPvariability$type <- factor(MSPvariability$type,levels = c("flow","price","econ"))

MSPvariability$CI95up <- as.numeric(MSPvariability$CI95up)*100
MSPvariability$CI95low <- as.numeric(MSPvariability$CI95low)*100
MSPvariability$mean <- MSPvariability$mean*100
ggplot(data=MSPvariability,aes(x=par,y=mean,fill=type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=CI95low,ymax=CI95up),color="black",position=position_dodge(.9),width=.4,size=2) +
  labs(x="Parameter", y="Variability (%)\n") +
  ggtitle(paste("MSP Variability (Flow)\n")) +
  scale_fill_Publication(labels=c("Process Parameters","Prices/Costs","Economic Assumptions")) +
  theme_Publication() +
  scale_x_discrete(labels=parNames) +
  theme(axis.title = element_text(size=25),
        axis.text.y = element_text(size=35),
        axis.text.x = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 45),
        legend.title = element_text(size=25,face="bold"),
        legend.text = element_text(size=20)) +
  guides(fill=guide_legend(
    title="Parameter Type",
    keywidth=0.3,
    keyheight=0.3,
    default.unit="inch"))
ggsave(file.path(output.dir,paste0("VariabilityMSP",".pdf")),width = 14, height = 8,dpi=500)
ggsave(file.path(output.dir,paste0("VariabilityMSP",".png")),width = 14, height = 8,dpi=500)

