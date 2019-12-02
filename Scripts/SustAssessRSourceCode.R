#SustAssessR source code specifying modules and functions for usage by the master implementation script. 
#This code requires process modules to have been seperately specified
#Author: Nick Matthews

#Define classes
#The process module class is the base class for "blocks" in the analysis
setClass("ProcessModule", slots=list(name="character", scaling="character", parameters="list",inputs="list",
                                 outputs="list",scaledInputs="list",scaledOutputs="list",
                                 foreground="list",background="list"))
#Used to calculate LCA results
setClass("LCAResult", slots=list(name="character",inputs="list",
                                 outputs="list",foreground="list",background="list",
                                 collated="list",collatedAverage="data.frame",aggregated="list",aggregatedAverage="data.frame",
                                 processFlows="data.frame",flowDists="data.frame",FermNumber="numeric"))
#Used to calculate costing results
setClass("CostingResult", slots=list(name="character",rawMats="list",rawMatsSum="numeric",
                                     utilities="list",utilitiesSum="numeric",waste="numeric",
                                     labour="numeric",otherFOC="numeric",capital="numeric",
                                     probComparison="numeric",combined="list",MSP="numeric",MSP_DF="data.frame",
                                     capCharge="numeric",breakdown="data.frame"))
#Used by allResultsProcess to bring all the results together and give relevant info
setClass("CombinedResults", slots=list(name="character",RawLCAresult="list",
                                       LCAresult="list",LCAresultsDF="data.frame",LCAHotspot="list",
                                       CostResult="list",MSP="list",CostBreakdown="list",probComparison="list",probComparisonDF="data.frame",
                                       CostBreakdownDF="data.frame",MSP_DF="data.frame",MSP_DFs="list",allResultDF="data.frame"))
#Used by nylonComparisons function for calculating comparisons in their nylon contexts
setClass("NylonResult", slots=list(collated="list",aggregated="data.frame"))

#####Utility functions#####
#Default impact methodology is ReCiPe2016, no other methodology is implemented yet
impactList <- function(name = 0, background=0,impactMeth = "ReCiPe2016" ) {
  #Set impact categories
  if(impactMeth == "ReCiPe2016") {
    impactCat <- c("Climate change, default, excl biogenic carbon", "Climate change, incl biogenic carbon", "Fine Particulate Matter Formation",
                   "Fossil depletion","Freshwater Consumption","Freshwater ecotoxicity","Freshwater Eutrophication","Human toxicity, cancer",
                   "Human toxicity, non-cancer","Ionizing Radiation","Land use","Marine ecotoxicity","Marine Eutrophication","Metal depletion",
                   "Photochemical Ozone Formation, Ecosystems","Photochemical Ozone Formation, Human Health","Stratospheric Ozone Depletion",
                   "Terrestrial Acidification","Terrestrial ecotoxicity")
  }
  #Initialise list
  temp <- list()
  #run list to create elements
  if(is.data.frame(background) & is.character(name)) {
    for(ii in 1:length(impactCat)) {
      temp[[rownames(background)[ii]]] <- background[rownames(background)==impactCat[ii],name]
    }
  } else {
    for(ii in 1:length(impactCat)) {
      temp[[impactCat[ii]]] <- 0
    }
  }
  #return output
  temp
}

#Function to scale a list by a given factor, for scaling background data
scaleList <- function(list,factor) {
  temp <- lapply(list,function(x) x*factor)
  temp
}

#Function to extract CF of certain kind - defaults to hierarchist perspective
getCF <- function(data.dir,name,impactCat,perspective="Hierarchist",CFs=NULL,compartment="none") {
  #First find or calculate characterisation factors
  if(is.null(CFs)) {
    recipe <- recipeCFs(data.dir)
  } else {
    recipe=CFs
  }
  #extract spreadsheet for specified impact cat
  temp <- recipe[[impactCat]]
  #Find names column
  temp$Name <- temp[,colnames(temp) %in% c("Name","name","Substance.name","Emitted.substance",
                                           "Name.in.ReCiPe","Mineral.resource","Fossil.resource")]
  #If the CF is specific to emmission compartment then that must be specified...
  if(sum(c("emission.compartment","Emission.compartment") %in% colnames(temp))==0) {
    return(temp[temp$Name==name,perspective])
  } else if(compartment=="none") {
    stop("Please specify emission compartment")
  } else {
    return(temp[temp$Name==name & temp[colnames(temp) %in% c("emission.compartment","Emission.compartment")]==compartment,perspective])
  }
}

#Function to process, collate and save ReCiPe files. If already exists, just load is
recipeCFs <- function(data.dir,recalculate=FALSE) {
  #Almost all times, the ReCiPeCFs file should already exist...
  if(file.exists(file.path(data.dir,"ReCiPeCFs.RData")) & recalculate == FALSE) {
    load(file.path(data.dir,"ReCiPeCFs.RData"))
  } else {
    #If it ReCiPeCFs file doesn't exist it needs to be created and saved
    files <- dir(file.path(data.dir,"CFs"))
    recipe <- list()
    for(ii in 1:length(files)) {
      temp <- read.csv(file.path(data.dir,"CFs",files[ii]),header = TRUE,skip=2,na.strings = c("","NA"),stringsAsFactors = FALSE)
      #Correct and ensure Hierachist perspective is there
      if("Hierarchist" %in% colnames(temp)) temp <- temp[!(is.na(temp$Hierarchist)|temp$Hierarchist=="#N/A"),]
      if("H" %in% colnames(temp)) {
        temp <- temp[!(is.na(temp$H)|temp$H=="#N/A"),]
        temp$Hierarchist <- temp$H
      }
      if("E" %in% colnames(temp)) {
        temp$Egalitarian <- temp$E
      }
      if("I" %in% colnames(temp)) {
        temp$Individualist <- temp$I
      }
      if("Hierarchic" %in% colnames(temp)) {
        temp <- temp[!(is.na(temp$Hierarchic)|temp$Hierarchic=="#N/A"),]
        temp$Hierarchist <- temp$Hierarchic
      }
      temp <- temp[,!is.na(temp[1,])]
      temp$Hierarchist <- as.numeric(gsub(",","",temp$Hierarchist))
      temp$Egalitarian <- as.numeric(gsub(",","",temp$Egalitarian))
      temp$Individualist <- as.numeric(gsub(",","",temp$Individualist))
      recipe[[gsub("ReCiPe.csv","",files[ii])]] <- temp
    }
    #Save file
    save(recipe,file=file.path(data.dir,"ReCiPeCFs.RData"))
  }
  #Either way, the relevant characterisation factors will returned
  return(recipe)
}

#Function to process general output data from GaBi
#At the moment it only works for recipe, can extract midpoint or endpoint indicators
cleanGabiResults <- function(raw,resultType="recipe",indicators = "midpoint",cleanIndicator = TRUE,removeQuant = TRUE) {
  #Select method based on result types
  if(resultType == "recipe") {
    #Remove top bits
    rawtable <- raw[(which(raw[,1] == "Quantities")-1):nrow(raw),]
    #remove all but last two columns
    rawtable <- rawtable[,max(grep("Quantities",rawtable[2,])):ncol(rawtable)]
    #Set colnames
    if(ncol(rawtable)==2) {
      rawtable[1,1:2] <- c("Indicator","Value")
    } else {
      rawtable[1,1] <- c("Indicator")
    }
    #Make column names first row
    colnames(rawtable) = as.character(unlist(rawtable[1, ]))
    rawtable = rawtable[-1, ] 
    #remove all zero values
    rawtable <- rawtable[rawtable[,2] != 0,]
    #Convert values column to numerics
    for(i in 2:ncol(rawtable)) {
      rawtable[,i] <- as.numeric(rawtable[,i])
    }
    #Extract only selected indicators
    if(indicators == "midpoint") {
      rawtable <- rawtable[grepl("Midpoint",rawtable$Indicator),]
    } else if(indicators == "endpoint") {
      rawtable <- rawtable[grepl("Endpoint",rawtable$Indicator),]
    }
    #Clean up indicators if requested
    if(cleanIndicator == TRUE) {
      rawtable$Indicator <- sapply(strsplit(rawtable$Indicator,"- "), function(x) x[2])
    }
    if(removeQuant == TRUE) {
      rawtable$Indicator <- sapply(strsplit(rawtable$Indicator," \\["), function(x) x[1])
    }
  } else {
    stop("No method for result type specified")
  }
  rawtable
}

#Function that takes a set of parameters and makes appropriate distributions
#Can cope with "none", "Unif", "Tri" and "Norm" at the moment
makeDists <- function(dists,samples=1000) {
  require(EnvStats,quietly = TRUE)
  #Initiate list
  parameters <- list()
  for(ii in 1:nrow(dists)) {
    if(dists$Dist[ii] == "none") {
      parameters[[rownames(dists)[ii]]] <- rep(dists$Value[ii],samples)
    } else if(dists$Dist[ii] == "Unif") {
      parameters[[rownames(dists)[ii]]] <- runif(samples,min=dists$Min[ii],max=dists$Max[ii])
    } else if(dists$Dist[ii] == "Tri") {
      if(dists$Min[ii] != dists$Value[ii] & dists$Max[ii] != dists$Value[ii]) {
        parameters[[rownames(dists)[ii]]] <- rtri(samples,min=dists$Min[ii],max=dists$Max[ii],mode=dists$Value[ii])
      } else if(dists$Min[ii] == dists$Max[ii]) {
        parameters[[rownames(dists)[ii]]] <- dists$Value[ii]
      } else if(dists$Min[ii] == dists$Value[ii]) {
        parameters[[rownames(dists)[ii]]] <- rtri(samples,min=dists$Min[ii],max=dists$Max[ii],mode=dists$Value[ii]+1E-10)
      } else if(dists$Max[ii] == dists$Value[ii]) {
        parameters[[rownames(dists)[ii]]] <- rtri(samples,min=dists$Min[ii],max=dists$Max[ii],mode=dists$Value[ii]-1E-10)
      }
    } else if(dists$Dist[ii] == "Norm") {
      parameters[[rownames(dists)[ii]]] <- rnorm(samples,mean=dists$Value[ii],sd=dists$SD[ii])
    }
  }
  tempDF <- parameters %>% Reduce(function(df1,df2) cbind(df1,df2),.)
  colnames(tempDF) <- names(parameters)
  return(tempDF)
}

#Function to extract all the key material and energy flows throughout the process
extractFlows <- function(tempRes) {
  require(magrittr, quietly=TRUE)
  #Flow to measure
  flowNames <- c("Steam","Electricity","MainFlowIn","MainFlowOut","AuxFlowOut","Sugar","Salt","Nitrogen","Water","Butanol","ButanolRec","NaOH")
  moduleNames <- names(tempRes@inputs)
  tempFlows <- list()
  for(ii in 1:length(moduleNames)) {
    #Set volume to minimum
    temp <- rep(0,length(flowNames))
    #Add in usage for every module
    for(jj in 1:length(temp)) {
      if(flowNames[jj] %in% names(tempRes@inputs[[ii]])) {
        temp[jj] <- temp[jj] + mean(tempRes@inputs[[ii]][[flowNames[jj]]])
      }
      if(flowNames[jj] %in% names(tempRes@outputs[[ii]])) {
        temp[jj] <- temp[jj] - mean(tempRes@outputs[[ii]][[flowNames[jj]]])
      }
    #Assign
    tempFlows[[moduleNames[ii]]] <- temp
    }
  }
  #
  tempDF <- tempFlows %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  colnames(tempDF) <- flowNames
  rownames(tempDF) <- moduleNames
  return(as.data.frame(tempDF))
}

#Function to take a list of distributions and summarise them
extractStats <- function(distributions) {
  require(EnvStats,quietly = TRUE)
  options(warn=-1) #Kurtosis function produces warning when run on fixed value (as for externally sourced values)
  tempDF <- data.frame(indicator=names(distributions),
                       mean=unlist(lapply(distributions,mean)),
                       min=unlist(lapply(distributions,min)),
                       max=unlist(lapply(distributions,max)),
                       CI95low=unlist(lapply(distributions,function(x) quantile(x,0.025))),
                       CI95up=unlist(lapply(distributions,function(x) quantile(x,0.975))),
                       var=unlist(lapply(distributions,var)),
                       sd=unlist(lapply(distributions,sd)),
                       variability=abs((unlist(lapply(distributions,max))-unlist(lapply(distributions,min)))/unlist(lapply(distributions,mean))),
                       skewness=unlist(lapply(distributions,sd)),
                       kurtosis=unlist(lapply(distributions,kurtosis)),
                       row.names=names(distributions),stringsAsFactors = FALSE)
  tempDF$kurtosis <- replace(tempDF$kurtosis,is.na(tempDF$kurtosis),0)
  options(warn=0)
  return(tempDF)
}

#Declare function for iteratively adding with mapply and returning in the desired format
mapplySpecial <- function(x,y,run){
  x[run]<-x[run]+y
  return(x)
}

#Function to fix all but one parateters for sensitivity analysis
fixParameters <- function(parameters,varPar,MCsims){
  #Create sample of right size (without replacement)
  tempPars <- parameters[sample(1:nrow(parameters),MCsims,replace=FALSE),]
  #Determine which row to fix (where to start)
  fixedRow <- sample(1:nrow(tempPars),1,replace=FALSE)
  #Fix all parameters but the tested one
  tempPars[,colnames(tempPars)!=varPar] <- do.call(rbind,replicate(MCsims,tempPars[fixedRow,colnames(tempPars)!=varPar],simplify = FALSE))
  return(tempPars)
}


#####Core functions#####
#All results integrated processing function
allResultsProcess <- function(data.dir,linkages,parameters,flowData,backgroundData,CFs,FU=1,targets,perspective,MCruns=1000) {
  #Create results object
  results <- new("CombinedResults")
  #Run for all perspectives and targets
  for(mm in 1:length(perspective)){
    for(kk in 1:length(targets)) {
      #Output calculation
      cat("Calculating all sustainability assessment results for",targets[kk],"with",nrow(parameters),"MC sims, using",perspective[mm],"Perspective\n")
      timings <- c()
      #Calculate LCA result
      start.time <- Sys.time()
      result <- LCAResultsProcess(data.dir,linkages,parameters,flowData,backgroundData,CFs,FU=FU,mol=targets[kk],perspective=perspective[mm],MCruns = MCruns)
      end.time <- Sys.time()
      timings["LCAProcessTime"] <- end.time - start.time
      #Extract LCA results
      results@RawLCAresult[[result@name]] <- result
      results@LCAresult[[result@name]] <- result@aggregatedAverage
      results@LCAHotspot[[result@name]] <- result@collatedAverage
      #Calculate costings results - only necessary for one perspective
      if(mm == 1){
        start.time <- Sys.time()
        results@CostResult[[result@name]] <- costCalc(data.dir,parameters,result,mol=targets[kk],MCruns=MCruns)
        end.time <- Sys.time()
        timings["CostingsTime"] <- end.time - start.time
        #Extract MSP and cost breakdown information
        results@MSP[[result@name]] <- results@CostResult[[result@name]]@MSP
        results@MSP_DFs[[result@name]] <- results@CostResult[[result@name]]@MSP_DF
        results@CostBreakdown[[result@name]] <- results@CostResult[[result@name]]@breakdown
      } else {
        timings["CostingsTime"] <- 0
      }
      cat("Timings were as follows:\n",paste(c("LCA:","Costings:"),c(round(timings,2))),"\n----\n")
    }
  }
  #Compile together LCA results
  results@LCAresultsDF <- results@LCAresult %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  results@CostBreakdownDF <- results@CostBreakdown %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  #Generate MSP dataframe
  results@MSP_DF <- results@MSP_DFs %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  
  #CombinedDF
  results@allResultDF <- rbind(results@LCAresultsDF,results@MSP_DF)
  
  return(results)
}

#This is the major script for calculating LCA results, using the modules specified seperately
LCAResultsProcess <- function(data.dir,linkages,parameters,flowData,backgroundData,CFs,FU=1,mol="Cadaverine",perspective="Hierarchist",MCruns=1000,
                              sensitivity=FALSE,variable="none") {
  #There are a number of packages required for this, they should already be loaded, but just in case
  require(magrittr, quietly=TRUE)

  #Establish new object for result
  tempRes <- new("LCAResult",name=paste0(perspective,mol))

  #Run for each MC iteration
  for(run in 1:nrow(parameters)) {
    #Load in background data based on perspective and geography
    background <- backgroundData[[paste0(parameters$geogs[run],"_Background_",perspective)]]
    #Choose linkages to use
    linkagesScen <- linkages[[paste0("linkages",parameters$WasteScenario[run])]]
    #Load in and extract flows - and create distributions
    #For US lignocellulosic there is a mix of glucose and xylose, so need to blend the data
    if(parameters$geogs[run] != "USLigno") {
      flows <- flowData[[paste0(mol,sugarScen[parameters$geogs[run]],parameters$DSP[run],"%DSPProcessed")]]
    } else {
      #If NREL scenario need to combine the flows proportional to different sugar contents
      flowsGlu <- flowData[[paste0(mol,"glucose",parameters$DSP[run],"%DSPProcessed")]]
      flowsXyl <- flowData[[paste0(mol,"xylose",parameters$DSP[run],"%DSPProcessed")]]
      flows <- cbind(flowsGlu[,c(1:4)],parameters$NRELGlucoseProp[run]*flowsGlu[,c(5:ncol(flowsGlu))] + parameters$NRELXyloseProp[run]*flowsXyl[,c(5:ncol(flowsXyl))])
    }

    #Now just need to select the flows, and correct a few of them...
    flowDists <- as.list(flows[,grepl(parameters$Switch[run],names(flows)) & grepl(parameters$OrgSpeed[run],names(flows)) & grepl(parameters$MaxYield[run],names(flows))])
    names(flowDists) <- flows[,1]
    #Ensure that all flow dists are correct - correct for excel based bugs
    flowDists <- flowDists %>% lapply(function(x) gsub(",","",x)) %>% lapply(as.numeric)
    #Add in flow dists for non-N in N fertiliser
    flowDists[["FertiliserMiscMassInput"]] <- flowDists$BioreactorNitrogenInput*parameters[["AmmoniumNitrogen"]][run]*parameters[[paste0("NonNitrogenRatio",parameters$NitrogenSource[run])]][run]
    #Apply NaOH correction
    flowDists[["pH_adjustmentNaOHInput"]] <- flowDists[["pH_adjustmentNaOHInput"]]*parameters$NaOHCorrection[run]
    #Add in butanol requirement
    flowDists[["ExtracterSolventInput"]] <- ((flowDists[["EvaporatorCakeOutput"]] + flowDists[["Distillation_step_2Heavy_OrganicsOutput"]] + 
                                                flowDists[["Distillation_step_2ProductOutput"]] + flowDists[["EvaporatorWaterOutput"]] + 
                                                flowDists[["EvaporatorWater_vaporOutput"]])/(1-parameters[["SolventReq"]][run])*parameters[["SolventReq"]][run])*parameters[["SolventLossRate"]][run]
    flowDists[["CakeWater"]] <- flowDists[["EvaporatorWaterOutput"]] + flowDists[["EvaporatorWater_vaporOutput"]]
    #Balance inputs and outputs
    flowDists$EvaporatorCakeOutput <- Reduce("+",flowDists[grepl("Input",names(flowDists))]) - 
      Reduce("+",flowDists[grepl("Output",names(flowDists)) & !grepl("Cake",names(flowDists))])
    #Set up flowDists data.frame in the result so we can record it
    if(run==1){
      tempRes@flowDists <- lapply(flowDists,function(x) rep(0,nrow(parameters))) %>% Reduce(function(df1,df2) cbind(df1,df2),.) %>% as.data.frame()
      colnames(tempRes@flowDists) <- names(flowDists)
    }
    tempRes@flowDists[run,] <- unlist(flowDists)
    #Run each linkage
    #Needs adaptation to iteratively append data
    for(ii in 1:nrow(linkagesScen)) {
      #If first one, do according to functional unit
      if(ii == 1) {
        temp <- do.call(linkagesScen[ii,"from"],args=list(data.dir,FU,run,parameters,flowDists,background,CFs,perspective,mol))
      } else {
        if(linkagesScen[ii,"type"]=="main") {
          temp <- do.call(linkagesScen[ii,"from"],args=list(data.dir,tempRes@inputs[[linkagesScen[ii,"to"]]][[linkagesScen[ii,"flow"]]][run],run,parameters,flowDists,background,CFs,perspective,mol))
        } else {
          temp <- do.call(linkagesScen[ii,"from"],args=list(data.dir,tempRes@outputs[[linkagesScen[ii,"to"]]][[linkagesScen[ii,"flow"]]][run],run,parameters,flowDists,background,CFs,perspective,mol))
        }
      }
      #If the module doesn't previously exist, need to create it
      if(is.null(tempRes@inputs[[linkagesScen[ii,"from"]]])) {
        tempRes@inputs[[linkagesScen[ii,"from"]]] <- lapply(temp@scaledInputs,function(x) rep(0,nrow(parameters)))
        tempRes@outputs[[linkagesScen[ii,"from"]]] <- lapply(temp@scaledOutputs,function(x) rep(0,nrow(parameters)))
        #Add foreground slot
        tempRes@foreground[[linkagesScen[ii,"from"]]] <- lapply(temp@foreground,function(x) rep(0,nrow(parameters)))
      }
      #Now add in the relevant results, The appropriate structures will now always exist, so just merge
      tempRes@inputs[[linkagesScen[ii,"from"]]] <- mapply(mapplySpecial,tempRes@inputs[[linkagesScen[ii,"from"]]],temp@scaledInputs,run,SIMPLIFY=FALSE)
      tempRes@outputs[[linkagesScen[ii,"from"]]] <- mapply(mapplySpecial,tempRes@outputs[[linkagesScen[ii,"from"]]],temp@scaledOutputs,run,SIMPLIFY=FALSE)
      #Add foreground results
      tempRes@foreground[[linkagesScen[ii,"from"]]] <- mapply(mapplySpecial,tempRes@foreground[[linkagesScen[ii,"from"]]],temp@foreground,run,SIMPLIFY=FALSE)

      #Add background results
      for(jj in 1:length(temp@background)) {
        #If background object doesn't previously exist, create it
        if(is.null(tempRes@background[[names(temp@background)[[jj]]]])) {
          tempRes@background[[names(temp@background)[[jj]]]] <- lapply(temp@background[[jj]],function(x) rep(0,nrow(parameters)))
        }
        #Add in data
        tempRes@background[[names(temp@background)[[jj]]]] <- mapply(mapplySpecial,tempRes@background[[names(temp@background)[[jj]]]],temp@background[[jj]],run,SIMPLIFY=FALSE)
        }
      }
    
  
    #If the scenario is USLigno need to run the lignocellulose step to account for impacts of sugar production
    if(parameters$geogs[run] %in% c("USLigno") & parameters$WasteScenario[run] != "NylonPolymerisation") {
      temp <- do.call(paste0("feedstock",parameters$geogs[run]),args=list(data.dir,tempRes@inputs$fermenter$Sugar[run],run,parameters,flowDists,background,CFs,perspective,mol))
      #record inputs and outputs
      #if doesn't previously exist, have to create it
      if(is.null(tempRes@inputs[[paste0("feedstock",parameters$geogs[run])]])) {
        tempRes@inputs[[paste0("feedstock",parameters$geogs[run])]] <- lapply(temp@scaledInputs,function(x) rep(0,nrow(parameters)))
        tempRes@outputs[[paste0("feedstock",parameters$geogs[run])]] <- lapply(temp@scaledOutputs,function(x) rep(0,nrow(parameters)))
        #Add foreground slot
        tempRes@foreground[[paste0("feedstock",parameters$geogs[run])]] <- lapply(temp@foreground,function(x) rep(0,nrow(parameters)))
      }
      #Add foreground results
      tempRes@inputs[[paste0("feedstock",parameters$geogs[run])]] <- mapply(mapplySpecial,tempRes@inputs[[paste0("feedstock",parameters$geogs[run])]],temp@scaledInputs,run,SIMPLIFY=FALSE)
      tempRes@outputs[[paste0("feedstock",parameters$geogs[run])]] <- mapply(mapplySpecial,tempRes@outputs[[paste0("feedstock",parameters$geogs[run])]],temp@scaledOutputs,run,SIMPLIFY=FALSE)
      #Add foreground results
      tempRes@foreground[[paste0("feedstock",parameters$geogs[run])]] <- mapply(mapplySpecial,tempRes@foreground[[paste0("feedstock",parameters$geogs[run])]],temp@foreground,run,SIMPLIFY=FALSE)
      
      #Add background results
      for(jj in 1:length(temp@background)) {
        #If doesn't previously exist, create it
        if(is.null(tempRes@background[[names(temp@background)[[jj]]]])) {
          tempRes@background[[names(temp@background)[[jj]]]] <- lapply(temp@background[[jj]],function(x) rep(0,nrow(parameters)))
        }
		#Add values for most recent dataset
        tempRes@background[[names(temp@background)[[jj]]]] <- mapply(mapplySpecial,tempRes@background[[names(temp@background)[[jj]]]],temp@background[[jj]],run,SIMPLIFY=FALSE)
      }
    }
    
    #Do not allow electricity to be in credit
    tempRes@background$Steam <- lapply(tempRes@background$Steam,function(x) {
       x[tempRes@background$Steam$`Climate change, default, excl biogenic carbon` < 0] <- 0
       return(x)})
  }
  #Collate results into one list of lists
  tempRes@collated <- c(tempRes@foreground,tempRes@background)
  
  
  #Average results into one dataframe
  tempRes@collatedAverage <- tempRes@collated %>%  lapply(lapply, mean) %>%
    simplify2array %>% data.frame %>% apply(2,as.numeric) %>% data.frame
  row.names(tempRes@collatedAverage) <- names(tempRes@collated[[1]])

  #Aggregate results into one list
  refNames <- names(tempRes@collated[[1]])
  for(ii in 1:length(tempRes@collated[[1]])) {
    temp <- lapply(tempRes@collated,function(x) x[[ii]])
    tempRes@aggregated[[refNames[ii]]] <- unlist(do.call(mapply,c(list(sum),temp)))
  }
  
  #Collate average aggregated results into one vector
  tempRes@aggregatedAverage <- extractStats(tempRes@aggregated)
  row.names(tempRes@aggregatedAverage) <- row.names(tempRes@collatedAverage)
  tempRes@aggregatedAverage$Scenario <- tempRes@name
  
  #Extract process flows
  if(unique(parameters$WasteScenario)[1] != "NylonPolymerisation"){
    tempRes@processFlows <- extractFlows(tempRes)
  }
  
  #Return output
  return(tempRes)
}

#Calculate costs
costCalc <- function(data.dir,parameters,LCAresult,mol="Cadaverine",MCruns=1000,sensitivity=FALSE,variable="none") {
  require(optiRum)

  #Create new costResult object
  costResult <- new("CostingResult",name=paste0(mol))
  #Establishe necessary vectors from start
  rawMats <- c("Sugar","Salt","Nitrogen","Water","Butanol","NaOH")
  utilities <- c("Steam","Electricity")
  #Set up empty vectors for recording
  for(rawMat in rawMats){
    costResult@rawMats[[rawMat]] <- rep(0,nrow(parameters))
  }
  for(utility in utilities){
    costResult@utilities[[utility]] <- rep(0,nrow(parameters))
  }
  costResult@waste <- rep(0,nrow(parameters))
  costResult@capital <- rep(0,nrow(parameters))
  costResult@labour <- rep(0,nrow(parameters))
  costResult@otherFOC <- rep(0,nrow(parameters))
  costResult@MSP <- rep(0,nrow(parameters))
  #Set up the discounted vectors
  costResult@combined$discountedCapital <- rep(0,nrow(parameters))
  costResult@combined$discountedInterest <- rep(0,nrow(parameters))
  costResult@combined$discountedLabour <- rep(0,nrow(parameters))
  costResult@combined$discountedOtherFOC <- rep(0,nrow(parameters))
  costResult@combined$discountedRawMats <- rep(0,nrow(parameters))
  costResult@combined$discountedUtilities <- rep(0,nrow(parameters))
  costResult@combined$discountedWaste <- rep(0,nrow(parameters))
  for(run in 1:nrow(parameters)) {
    ##Raw Materials##
    #Raw materials - at a per kg basis
    for(rawMat in rawMats) {
      #Set volume to minimum
      totalVol <- 0
      #Add in usage for every module
      for(ii in which(sapply(LCAresult@inputs,function(x) rawMat %in% names(x)))) {
        totalVol <- totalVol + LCAresult@inputs[[ii]][[rawMat]][run]
      }
      #Extract costs to object - correction if it's lignocellulosic sugar
      if(rawMat == "Sugar" & parameters$geogs[run] %in% c("USLigno")) {
        costResult@rawMats[[rawMat]][run] <- totalVol*parameters$LignoSugar[run]*1E8
      } else {
        costResult@rawMats[[rawMat]][run] <- totalVol*parameters[[rawMat]][run]*1E8
      }
    }
    #Correct nitrogen calculation seperately
    costResult@rawMats$Nitrogen[run] <- costResult@rawMats$Nitrogen[run]*(1+parameters[[paste0("NonNitrogenRatio",parameters$NitrogenSource[run])]][run])
    #Aggregate costs together
    costResult@rawMatsSum[run] <- sum(unlist(lapply(costResult@rawMats,"[[",run)))
    
    ##Utilties##
    #Do similar for utilities
    for(utility in utilities) {
      #Set volume to minimum
      totalVol <- 0
      #Add in usage for every module
      for(ii in which(sapply(LCAresult@inputs,function(x) utility %in% names(x)))) {
        totalVol <- totalVol + LCAresult@inputs[[ii]][[utility]][run]
      }
      for(ii in which(sapply(LCAresult@outputs,function(x) utility %in% names(x)))) {
        totalVol <- totalVol - LCAresult@outputs[[ii]][[utility]][run]
      }
      #if steam correct the amount (MJ to kg) and don't allow to be in credit
      if(utility == "Steam") {
        #MJ -> kg
        totalVol <- totalVol/2.1
        #Do not allow steam to be in credit as not useable
        totalVol[totalVol < 0] <- 0
      }
      #Extract costs to object
      costResult@utilities[[utility]][run] <- totalVol*parameters[[utility]][run]*1E8
    }
    costResult@utilitiesSum[run] <- sum(unlist(lapply(costResult@utilities,"[[",run)))
  
    ##Waste treatment##
    #Set bolume to minimum
    totalVol <- 0
    #Add in usage for every module
    for(ii in which(sapply(LCAresult@outputs,function(x) "Waste" %in% names(x)))) {
      totalVol <- totalVol + LCAresult@outputs[[ii]][["Waste"]][run]
    }
    costResult@waste[run] <- totalVol * parameters[["Waste"]][run] * 1E8
  
    ##Capital costs##
    costResult@capital[run] <- parameters$BioethanolPlant[run]*(LCAresult@flowDists$GeneralNumber_of_fermentersParameter[run]/
                                                                  parameters$BaseFermSize[run])^parameters$CapitalScaling[run]*1E6
    
    #####Operating costs#####
    totLab <- parameters$Labour[run]*(LCAresult@flowDists$GeneralNumber_of_fermentersParameter[run]/
                                            parameters$BaseFermSize[run])^parameters$LabourScaling[run]
    TCI <- parameters$BioethanolPlant[run]*(LCAresult@flowDists$GeneralNumber_of_fermentersParameter[run]/
                                              parameters$BaseFermSize[run])^parameters$CapitalScaling[run]
    #Laborotaory consideration
    tax <- TCI * parameters$Tax[run]
    maintenance <- TCI * parameters$Maintenance[run]
    FOC <- tax + totLab + maintenance
    RDMarketing <- FOC*parameters$RDMarketing[run]
    overheads <- totLab*parameters$overheads[run]
    costResult@labour[run] <- totLab * 1E6
    costResult@otherFOC[run] <- (maintenance + tax + RDMarketing + overheads) * 1E6
  
    ##Now calculate MSP##
    years <- seq(1,parameters$times[run],1)
    combinedCosts <- list()
    discountFactor <- 1/(1+parameters$discountRate[run])^years
    annualCost <- costResult@rawMatsSum[run] + costResult@utilitiesSum[run] + costResult@waste[run] + costResult@labour[run] + costResult@otherFOC[run]
    
    combinedCosts[["discountedRawMats"]] <- costResult@rawMatsSum[run] * discountFactor
    combinedCosts$discountedRawMats[1] <- combinedCosts$discountedRawMats[1]*0.5 #First year on 50% production
    combinedCosts[["discountedUtilities"]] <- costResult@utilitiesSum[run]*discountFactor
    combinedCosts$discountedUtilities[1] <- combinedCosts$discountedUtilities[1]*0.5 #First year on 50% production
    combinedCosts[["discountedWaste"]] <- costResult@waste[run] *discountFactor
    combinedCosts$discountedWaste[1] <- combinedCosts$discountedWaste[1]*0.5 #First year on 50% production
    combinedCosts[["discountedLabour"]] <- costResult@labour[run] *discountFactor
    combinedCosts$discountedLabour[1] <- combinedCosts$discountedLabour[1]*0.75 #First year on 75% labour
    combinedCosts[["discountedOtherFOC"]] <- costResult@otherFOC[run] *discountFactor
    combinedCosts[["discountedAnnualCost"]] <- combinedCosts$discountedRawMats+combinedCosts$discountedUtilities+
                                                            combinedCosts$discountedWaste+combinedCosts$discountedLabour+
                                                            combinedCosts$discountedOtherFOC
    #Calculate loan amount
    capitalTotal <- TCI*1E6
    depreciation <- c(0.14,	0.2449,	0.1749,	0.1249,	0.0893,	0.0893,	0.0892,	0.0475)
    interestRate <- parameters$interestRate[run]
    loanAcru <- (((capitalTotal*0.08*0.6)*(1+interestRate) + capitalTotal*0.6*0.6)*(1+interestRate) + capitalTotal*0.32*0.6)*(1+interestRate)
    #Calculate repayments required
    repayments <- PMT(interestRate,round(parameters$RepaymentPeriod[run]),loanAcru)
    interestPayment <- c()
    loanRemaining <- loanAcru
    for(i in 1:round(parameters$RepaymentPeriod[run])) {
      interestPayment[i] <- loanRemaining*interestRate
      loanRemaining <- loanRemaining + repayments + interestPayment[i]
    }
    interestPayment <- c(interestPayment,rep(0,length(years)-round(parameters$RepaymentPeriod[run])))
    combinedCosts[["discountedInterest"]] <- discountFactor*interestPayment
    capital <- c(depreciation*capitalTotal,rep(0,length(years)-length(depreciation)))
    combinedCosts[["discountedCapital"]] <- discountFactor*capital
    combinedCosts[["annualisedCosts"]] <- combinedCosts$discountedAnnualCost+combinedCosts$discountedCapital+combinedCosts$discountedInterest
    for(MSPtemp in seq(0.1,20,0.05)) {
      income <- rep(MSPtemp*1E8,parameters$times[run])
      #Adjust for reduced time in year 1
      income[1] <- MSPtemp*0.5*1E8
      #Subtract costs from income
      incomeDiscounted <- income*discountFactor
      discountedGrossIncome <- incomeDiscounted-combinedCosts$annualisedCosts
      #Apply income tax rate if making money
      discountedNetIncome <- unlist(lapply(discountedGrossIncome,function(x) unlist(mapply(function(y,z){if(y>0) y*(1-z) else y},x,parameters$incomeTax[run]))))
      totalNetIncome <- sum(discountedNetIncome)
      #If total netincome more than 0, then this is the MSP
      if(totalNetIncome > 0){
        costResult@MSP[run] <- MSPtemp
        break()
      }
    }
    
    #Record the cost results
    costResult@combined$discountedCapital[run] <- sum(combinedCosts$discountedCapital)
    costResult@combined$discountedInterest[run] <- sum(combinedCosts$discountedInterest)
    costResult@combined$discountedLabour[run] <- sum(combinedCosts$discountedLabour)
    costResult@combined$discountedOtherFOC[run] <- sum(combinedCosts$discountedOtherFOC)
    costResult@combined$discountedRawMats[run] <- sum(combinedCosts$discountedRawMats)
    costResult@combined$discountedUtilities[run] <- sum(combinedCosts$discountedUtilities)
    costResult@combined$discountedWaste[run] <- sum(combinedCosts$discountedWaste)
  }
  #Create breakdown table
  costResult@breakdown <- extractStats(costResult@combined)
  costResult@breakdown$Area <- c("CapitalDepreciation","InterestPayments","Labour","OtherFOC","RawMaterials","Utilities","Waste")
  costResult@breakdown$Scenario <- LCAresult@name
  #Extract summary of stats for MSP
  costResult@MSP_DF <- extractStats(list(MSP=costResult@MSP))
  costResult@MSP_DF$Scenario <- costResult@name
  #Return object
  return(costResult)
}

#Function to do normalisation
normalisation <- function(data.dir,result,perspectives=c("Egalitarian","Hierarchist","Individualist"),targets=c("Cadaverine","Putrescine"),resultsUsed="none",type="none") {
  #First load in recipe normalisations
  norm <- read.csv(file.path(data.dir,"RecipeNormalisation250119_2.csv"),header = TRUE,stringsAsFactors = FALSE)
  #Extract all the scenarios to look at and extract necessary information - if this is one of our results
  if(class(result)[1] == "CombinedResults") {
    allScen <- names(result@RawLCAresult)
    resultsUsed <- which(grepl(paste(targets,collapse="|"),allScen))
    resDF <- result@LCAresultsDF[grepl(paste(targets,collapse="|"),result@LCAresultsDF$Scenario),]
    #Extract perspective and geographical info
    resDF$Perspective <- "none"
    resDF$Perspective[grepl(paste("Hierarchist",collapse="|"),resDF$Scenario)] <- "Hierarchist"
    resDF$Perspective[grepl(paste("Egalitarian",collapse="|"),resDF$Scenario)] <- "Egalitarian"
    resDF$Perspective[grepl(paste("Individualist",collapse="|"),resDF$Scenario)] <- "Individualist"
  } else if(type=="nylon") {
    #If nylon comparison results, process this way
    allScen <- unique(result$key)
    resultsUsed <- 1:length(allScen)
    resDF <- result
    row.names(resDF) <- NULL
  } else if(type == "comparison"){
    #If not a self-calculated result, require a different kind of processing
    allScen <- unique(result$Scenario)
    resultsUsed <- 1:length(allScen)
    resDF <- result
    resDF$result <- resDF$Value
    resDF$CIup <- resDF$Value
    resDF$CIlow <- resDF$Value
    resDF$indicator <- resDF$Indicator
  } else {
    stop("Unrecognised data type")
  }
  
  #Establish dataframes
  tempList <- list()
  for(ii in resultsUsed) {
    #Extract necessary base results
    if(type=="nylon"){
      baseDF <- resDF[resDF$key==allScen[ii],]
    } else {
      baseDF <- resDF[resDF$Scenario==allScen[ii],]
    }
    #Select the necessary normalisation factors
    tempDF <- norm[norm$Country=="Glo",c("Name","Country","Midpoint","Endpoint","unit",unique(baseDF$Perspective))]
    #Get the relevant normalisation factors
    tempDF$NormScore <- tempDF[,unique(baseDF$Perspective)]
    tempDF[,unique(baseDF$Perspective)] <- NULL
    #extract the necessary background info
    tempDF$Perspective <- unique(baseDF$Perspective)
    if(type=="nylon") {
      tempDF$nylon <- unique(baseDF$nylon)
      tempDF$key <- unique(baseDF$key)
    }
    #Do the normalisation - matching to find the info
    tempDF$EndScoreMean <- baseDF$mean[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    tempDF$EndScoreMax <- baseDF$min[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    tempDF$EndScoreMin <- baseDF$max[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    tempDF$EndScoreVariability <- baseDF$variability[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    tempDF$EndScoreCI95up <- baseDF$CI95up[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    tempDF$EndScoreCI95low <- baseDF$CI95low[match(tempDF$Midpoint,baseDF$indicator)]*tempDF$NormScore
    if(type=="nylon"){
      tempList[[unique(baseDF$key)]] <- tempDF
    } else {
      tempList[[unique(baseDF$Scenario)]] <- tempDF
    }
  }
  resultsDF <- tempList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
  
  return(resultsDF)
}

#Function which extracts particular comparison and merges it together
nylonComparisons <- function(result,comparisonFile,parameters,MolPar,linkages,flowData,backgroundData,CFs,data.dir,perspective="Hierarchist",mainNylon="nylon46",refNylon="Nylon66",main="Putrescine",reactant="AdipicAcid_RoW",reference=c("AdipicAcid_RoW","HMDA_Glo")) {
  cat("Calculating nylon comparisons for",mainNylon,"vs",refNylon,"with",perspective,"Perspective\n")
  start.time=Sys.time()
  require(stringr, quietly=TRUE)
  #Create nylon result data
  nylonRes <- new("NylonResult")
  #Select scenario datasets to use
  allScen <- names(result@RawLCAresult)
  resultsUsed <- which(grepl(paste(main,collapse="|"),allScen) & 
                        grepl(paste(perspective,collapse="|"),allScen))
  #Establish dataframe for stoichiometry calculations
  tempDF <- data.frame(mol=c(main,reactant,reference),
                       key=c(rep("main",length(main)),rep("main",length(reactant)),rep("reference",length(reference))),
                       stringsAsFactors = FALSE)
  #Find Mr based on their names
  tempDF$Mr <- rep(0,nrow(tempDF))
  for(ii in 1:nrow(tempDF)) {
    tempDF$Mr[ii] <- MolPar[paste0(str_split(tempDF$mol[ii],"_")[[1]][1],"Mr"),"Value"]
  }
  #Calculate ratios from Mr, assuming 1:1 ratio of nylons
  tempDF$ratio <- rep(0,nrow(tempDF))
  for(component in unique(tempDF$key)) {
    total <- sum(tempDF$Mr[tempDF$key==component])
    tempDF$ratio[tempDF$key==component] <- tempDF$Mr[tempDF$key==component]/total
  }
  #Extract impact lists for each of the components
  tempList <- list()
  #Add in core results
  for(ii in resultsUsed) {
    tempList[[allScen[ii]]] <- scaleList(result@RawLCAresult[[allScen[ii]]]@aggregated,tempDF$ratio[1])
  }
  
  #Now add in comparisons
  for(ii in 2:nrow(tempDF)){
    if(tempDF$mol[ii] == "SebacicAcid_average"){
      tempListItem <- impactList()
      for(impactCat in names(tempListItem)){
        tempListItem[[impactCat]] <- as.numeric(comparisonFile[impactCat,parameters$SebacicAcidComp])
      }
      tempList[[paste0(tempDF$mol[ii],"_",tempDF$key[ii])]] <- scaleList(tempListItem,tempDF$ratio[ii])
    } else if(tempDF$mol[ii] == "AdipicAcid_average"){
      tempListItem <- impactList()
      for(impactCat in names(tempListItem)){
        tempListItem[[impactCat]] <- as.numeric(comparisonFile[impactCat,parameters$AdipicAcidComp])
      }
      tempList[[paste0(tempDF$mol[ii],"_",tempDF$key[ii])]] <- scaleList(tempListItem,tempDF$ratio[ii])
    } else {
      tempList[[paste0(tempDF$mol[ii],"_",tempDF$key[ii])]] <- scaleList(impactList(paste0(tempDF$mol[ii]),comparisonFile),tempDF$ratio[ii])
    }
  }
  
  #
  parametersMod <- parameters
  parametersMod$WasteScenario <- "NylonPolymerisation"
  
  #Compute processing results - use LCA results process
  polyResRaw <- LCAResultsProcess(data.dir,linkages,parametersMod,flowData,backgroundData,CFs,FU=1,perspective=perspective,MCruns=nrow(parametersMod))
  polyAggregated <- polyResRaw@aggregatedAverage
  polyAggregated$mol <- "Polymerisation"
  polyAggregated$Scenario <- NULL
  polyAggregatedMain <- polyAggregated
  polyAggregatedMain$key <- "main"
  polyAggregatedReference <- polyAggregated
  polyAggregatedReference$key <- "reference"
  
  #One collated uncertain dataset
  #Need to add all the distributions together
  refNames <- names(tempList[[1]])
  for(ii in 1:length(tempList[[1]])) {
    temp <- lapply(tempList[1:2],function(x) x[[ii]])
    temp$polymerisation <- polyResRaw@aggregated[[refNames[ii]]]
    nylonRes@collated[[mainNylon]][[refNames[ii]]] <- unlist(do.call(mapply,c(list(sum),temp)))
  }
  for(ii in 1:length(tempList[[1]])) {
    temp <- lapply(tempList[3:4],function(x) x[[ii]])
    temp$polymerisation <- polyResRaw@aggregated[[refNames[ii]]]
    nylonRes@collated[[refNylon]][[refNames[ii]]] <- unlist(do.call(mapply,c(list(sum),temp)))
  }
  
  #Extract stats
  mon1Main <- extractStats(tempList[[allScen[resultsUsed]]])
  mon1Main$mol <- tempDF$mol[1]
  mon1Main$key <- tempDF$key[1]
  mon2Main <- extractStats(tempList[[paste0(tempDF$mol[tempDF$key=="main"][2],"_main")]])
  mon2Main$mol <- tempDF$mol[2]
  mon2Main$key <- tempDF$key[2]

  mon1Reference <- extractStats(tempList[[paste0(tempDF$mol[tempDF$key=="reference"][1],"_reference")]])
  mon1Reference$mol <- tempDF$mol[3]
  mon1Reference$key <- tempDF$key[3]

  mon2Reference <- extractStats(tempList[[paste0(tempDF$mol[tempDF$key=="reference"][2],"_reference")]])
  mon2Reference$mol <- tempDF$mol[4]
  mon2Reference$key <- tempDF$key[4]
  
  #Need to sum up each one
  mainCollated <- mon1Main
  mainCollated$mol <- "Collated"
  mainCollated[,2:11] <- 0
  mainCollated$mean <- mon1Main$mean + mon2Main$mean + polyAggregatedMain$mean
  mainCollated$min <- mon1Main$min + mon2Main$min + polyAggregatedMain$min
  mainCollated$max <- mon1Main$max + mon2Main$max + polyAggregatedMain$max
  mainCollated$CI95up <- mon1Main$CI95up + mon2Main$CI95up + polyAggregatedMain$CI95up
  mainCollated$CI95low <- mon1Main$CI95low + mon2Main$CI95low + polyAggregatedMain$CI95low
  mainCollated$variability <- (mainCollated$max-mainCollated$min)/mainCollated$mean
  
  referenceCollated <- mon1Reference
  referenceCollated$mol <- "Collated"
  referenceCollated[,2:11] <- 0
  referenceCollated$mean <- mon1Reference$mean + mon2Reference$mean + polyAggregatedReference$mean
  referenceCollated$min <- mon1Reference$min + mon2Reference$min + polyAggregatedReference$min
  referenceCollated$max <- mon1Reference$max + mon2Reference$max + polyAggregatedReference$max
  referenceCollated$CI95up <- mon1Reference$CI95up + mon2Reference$CI95up + polyAggregatedReference$CI95up
  referenceCollated$CI95low <- mon1Reference$CI95low + mon2Reference$CI95low + polyAggregatedReference$CI95low
  referenceCollated$variability <- (referenceCollated$max-referenceCollated$min)/referenceCollated$mean
  
  nylonRes@aggregated <- rbind(mon1Main,mon2Main,polyAggregatedMain,mainCollated,mon1Reference,mon2Reference,polyAggregatedReference,referenceCollated)
  nylonRes@aggregated$nylon[nylonRes@aggregated$key=="main"] <- mainNylon
  nylonRes@aggregated$nylon[nylonRes@aggregated$key=="reference"] <- refNylon
  nylonRes@aggregated$Perspective <- perspective
  nylonRes@aggregated$Scenario <- paste0(nylonRes@aggregated$nylon,nylonRes@aggregated$Perspective)
  
  end.time=Sys.time()
  cat("Time taken:",end.time-start.time,"\n------\n")
  #Return final result
  return(nylonRes)
  
}

#Run sensitivity test
sensTest <- function(data.dir,baseResult=NULL,sensTestPars,linkages,parameters,CFs,target="Cadaverine",perspective="Hierarchist",type="MinMax",MCruns=1000,SensSims=10,MultiStart=1) {
  cat("Calculating variability using",type,"method for",target,"with",perspective,"perspective.\n","Type:",type,"MultiStarts:",MultiStart,"\n-----\n")

  #First identify overall variability, only recalculate results if necessary
  if(is.null(baseResult)){
    result <- LCAResultsProcess(data.dir,linkages,parameters,flowData,backgroundData,CFs,FU=1,mol=target,perspective=perspective,MCruns = MCruns)
    costResult <- costCalc(data.dir,parameters,result,mol=target,MCruns=MCruns)
  } else {
    result <- baseResult@RawLCAresult[[paste0(perspective,target)]]
    costResult <- baseResult@CostResult[[grep(target,names(baseResult@CostResult))[1]]]
  }
  
  #Create sensitivity dataframe
  sensitivity <- rbind(result@aggregatedAverage,costResult@MSP_DF)
  sensitivity$indicator <- as.character(sensitivity$indicator)
  sensitivity$par <- "overall"
  sensitivity$type <- "overall"
  sensitivity$mean <- sensitivity$variability
  sensitivity$min <- sensitivity$variability
  sensitivity$max <- sensitivity$variability
  sensitivity[,5:11] <- 0
  indicators <- sensitivity$indicator 
  #Run sensitivity - just single start
  for(ii in 1:nrow(sensTestPars)) {
    start.time <- Sys.time()
    cat("Calculating sensitivity for parameter",sensTestPars$name[ii],"\n")
    varList <- list()
    #Identify ranges of values for parameter
    if(type=="MinMax"){
      if(is.character(parameters[,sensTestPars$name[ii]])){
        #If character fix to unique values
        parValues <- unique(parameters[,sensTestPars$name[ii]])
      } else if(is.numeric(parameters[,sensTestPars$name[ii]])){
        minPar <- min(parameters[sensTestPars$name[ii]])
        maxPar <- max(parameters[sensTestPars$name[ii]])
        parValues <- c(minPar,maxPar)
        #Need ligno parvalues
        if(sensTestPars$name[ii] == "Sugar"){
          minPar <- min(parameters["LignoSugar"])
          maxPar <- max(parameters["LignoSugar"])
          parValuesLigno <- c(minPar,maxPar)
        }
      }
    }
    for(jj in 1:MultiStart){
      #Extract the necessary parameter set
      if(type=="MinMax") {
        tempParameters <- fixParameters(parameters,varPar=sensTestPars$name[ii],MCsims=length(parValues))
        tempParameters[sensTestPars$name[ii]] <- parValues
        if(sensTestPars$name[ii] == "Sugar"){tempParameters["LignoSugar"] <- parValues}
      } else {
        tempParameters <- fixParameters(parameters,varPar=sensTestPars$name[ii],MCsims=SensSims)
      }
      
      #Run LCA and costing calculation
      tempResult <- LCAResultsProcess(data.dir,linkages,tempParameters,flowData,backgroundData,CFs,FU=1,mol=target,perspective=perspective,
                                      MCruns = nrow(tempParameters))
      costResult <- costCalc(data.dir,tempParameters,tempResult,mol=target,MCruns=nrow(tempParameters))
      tempVar <- rbind(tempResult@aggregatedAverage,costResult@MSP_DF)
      tempVar$indicator <- as.character(tempVar$indicator)
      tempVar$par <- sensTestPars$name[ii]
      varList[[jj]] <- tempVar
    }
    allVar <- varList %>% Reduce(function(df1,df2) rbind(df1,df2),.)
    if(sum(allVar$variability)>0) {
      #Collate all variability into set of lists
      allVarList <- lapply(unique(allVar$indicator),function(x) allVar$variability[allVar$indicator==x])
      names(allVarList) <- unique(allVar$indicator)
      tempVar <- extractStats(allVarList)
      tempVar$Scenario <- unique(allVar$Scenario)[1]
      tempVar$par <- unique(allVar$par)[1]
      tempVar$type <- sensTestPars$type[ii]
      sensitivity <- rbind(sensitivity,tempVar)
    } else {
      cat("Parameter",sensTestPars$name[ii],"has no sensitivity\n")
    }
    warnings()
    end.time <- Sys.time()
    cat("Time taken:",end.time-start.time,"\n------\n")
  }
  #Now return the result
  return(sensitivity)
}


