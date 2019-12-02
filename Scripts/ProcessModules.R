#Source code containing module specifications for process model
#Author: Nick Matthews

#Define polymerisation module
polymerisation <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="polymerisation",scaling="Plastic",foreground=impactList())
  #Define parameters
  module@parameters[["Transport"]] <- parameters[["PlasticTransport"]][run]
  module@parameters[["ElecRequired"]] <- parameters[["PlasticElectricity"]][run]/3.6 #Must convert to kwh
  module@parameters[["SteamRequired"]] <- parameters[["PlasticThermal"]][run]
  module@parameters[["PlasticProduced"]] <- 1

  #Define inputs
  module@inputs[["Electricity"]] <- module@parameters[["ElecRequired"]]
  module@inputs[["Steam"]] <- module@parameters[["SteamRequired"]]
  module@inputs[["Transport"]] <- module@parameters[["Transport"]]
  
  #Define outputs
  module@outputs[["Plastic"]] <- module@parameters[["PlasticProduced"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Plastic"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Electricity"]] <- scaleList(impactList(paste0("Electricity",parameters$energySource[run]),background),module@scaledInputs[["Electricity"]])
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])
  module@background[["Transport"]] <- scaleList(impactList("Transport",background),module@scaledInputs[["Transport"]])
  
  #Return object
  return(module)
  
  
}

#Define lignocellulose scenario processing impacts
feedstockUSLigno <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Ligno Sugar Production",scaling="Sugar",foreground=impactList())
  #Define parameters
  module@parameters[["CornStover"]] <- parameters[["NRELCornStover"]][run]
  module@parameters[["SulfuricAcid"]] <- parameters[["NRELSulfuricAcid"]][run]
  module@parameters[["Lime"]] <- parameters[["NRELLime"]][run]
  module@parameters[["Water"]] <- parameters[["NRELWater"]][run]
  module@parameters[["ElecGen"]] <- parameters[["NRELElecGen"]][run]
  module@parameters[["Waste"]] <- parameters[["NRELAsh"]][run]
  module@parameters[["SugarGen"]] <- parameters[["NRELSugarGen"]][run]
  module@parameters[["CornSteepLiqour"]] <- parameters[["NRELCornSteepLiqour"]][run]
  module@parameters[["NaOH"]] <- parameters[["NRELNaOH"]][run]
  module@parameters[["CornSteepLiqour"]] <- parameters[["NRELCornSteepLiqour"]][run]
  module@parameters[["Ammonia"]] <- parameters[["NRELAmmonia"]][run]
  module@parameters[["Glucose"]] <- parameters[["NRELGlucose"]][run]
  module@parameters[["AmmoniumSulfate"]] <- parameters[["NRELHostNutrients"]][run]
  module@parameters[["SoybeanOil"]] <- parameters[["NRELCornOil"]][run]
  module@parameters[["SO2"]] <- parameters[["NRELSO2"]][run]
  module@parameters[["SO2Emit"]] <- parameters[["NRELEmitSO2"]][run]
  module@parameters[["CO2Emit"]] <- parameters[["NRELEmitCO2"]][run]
  module@parameters[["COEmit"]] <- parameters[["NRELEmitCO"]][run]
  module@parameters[["NO2Emit"]] <- parameters[["NRELEmitNO2"]][run]
  module@parameters[["CH4Emit"]] <- parameters[["NRELEmitCH4"]][run]
  
  #Define inputs
  module@inputs[["Lime"]] <- module@parameters[["Lime"]]
  module@inputs[["CornStover"]] <- module@parameters[["CornStover"]]
  module@inputs[["Water"]] <- module@parameters[["Water"]]
  module@inputs[["SulfuricAcid"]] <- module@parameters[["SulfuricAcid"]]
  module@inputs[["CornSteepLiqour"]] <- module@parameters[["CornSteepLiqour"]]
  module@inputs[["NaOH"]] <- module@parameters[["NaOH"]]
  module@inputs[["Ammonia"]] <- module@parameters[["Ammonia"]]
  module@inputs[["Glucose"]] <- module@parameters[["Glucose"]]
  module@inputs[["AmmoniumSulfate"]] <- module@parameters[["AmmoniumSulfate"]]
  module@inputs[["SoybeanOil"]] <- module@parameters[["SoybeanOil"]]
  module@inputs[["SO2"]] <- module@parameters[["SO2"]]
  
  #Define outputs
  module@outputs[["Electricity"]] <- module@parameters[["ElecGen"]]
  module@outputs[["Sugar"]] <- module@parameters[["SugarGen"]]
  module@outputs[["Waste"]] <- module@parameters[["Waste"]]
  module@outputs[["CO2Emit"]] <- module@parameters[["CO2Emit"]]
  module@outputs[["COEmit"]] <- module@parameters[["COEmit"]]
  module@outputs[["SO2Emit"]] <- module@parameters[["SO2Emit"]]
  module@outputs[["NO2Emit"]] <- module@parameters[["NO2Emit"]]
  module@outputs[["CH4Emit"]] <- module@parameters[["CH4Emit"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Sugar"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Electricity"]] <- scaleList(impactList(paste0("Electricity",parameters$energySource[run]),background),-module@scaledOutputs[["Electricity"]])
  module@background[["Lime"]] <- scaleList(impactList("Lime",background),module@scaledInputs[["Lime"]])
  module@background[["CornStover"]] <- scaleList(impactList("CornStover",background),module@scaledInputs[["CornStover"]])
  module@background[["Glucose"]] <- scaleList(impactList("Sugar",background),module@scaledInputs[["Glucose"]])
  module@background[["Water"]] <- scaleList(impactList("Water",background),module@scaledInputs[["Water"]])
  module@background[["SulfuricAcid"]] <- scaleList(impactList("SulfuricAcid",background),module@scaledInputs[["SulfuricAcid"]])
  module@background[["NaOH"]] <- scaleList(impactList("SodiumHydroxide",background),module@scaledInputs[["NaOH"]])
  module@background[["CornSteepLiqour"]] <- scaleList(impactList("CornSteepLiqour",background),module@scaledInputs[["CornSteepLiqour"]])
  module@background[["SoybeanOil"]] <- scaleList(impactList("SoybeanOil",background),module@scaledInputs[["SoybeanOil"]])
  module@background[["SO2"]] <- scaleList(impactList("SO2",background),module@scaledInputs[["SO2"]])
  module@background[["Ammonia"]] <- scaleList(impactList("Ammonia",background),module@scaledInputs[["Ammonia"]])
  module@background[["AmmoniumSulfate"]] <- scaleList(impactList("AmmoniumSulfate",background),module@scaledInputs[["AmmoniumSulfate"]])
  module@background[["Waste"]] <- scaleList(impactList("WasteTreatment2",background),module@scaledOutputs[["Waste"]])
  
  
  #Define foreground impacts
  #For CO2
  module@foreground[["Climate change, incl biogenic carbon"]] <- module@scaledOutputs[["CO2Emit"]]*getCF(data.dir,"Carbon dioxide","GlobalWarming",perspective,CFs)
  #For CH4
  module@foreground[["Climate change, incl biogenic carbon"]] <- module@scaledOutputs[["CH4Emit"]]*getCF(data.dir,"Methane","GlobalWarming",perspective,CFs)
  module@foreground[["Climate change, default, excl biogenic carbon"]] <- module@scaledOutputs[["CH4Emit"]]*getCF(data.dir,"Methane","GlobalWarming",perspective,CFs)
  #For SO2
  module@foreground[["Fine Particulate Matter Formation"]] <- module@scaledOutputs[["SO2Emit"]]*getCF(data.dir,"SO2","ParticulateMatter",perspective,CFs)
  module@foreground[["Terrestrial Acidification"]] <- module@scaledOutputs[["SO2Emit"]]*getCF(data.dir,"SO2","TerrestrialAcidification",perspective,CFs)
  #For NO2
  module@foreground[["Photochemical Ozone Formation, Human Health"]] <- module@scaledOutputs[["NO2Emit"]]*getCF(data.dir,"nitrogen oxides (as nitrogen dioxide)","OzoneFormationHuman",perspective,CFs)
  module@foreground[["Photochemical Ozone Formation, Ecosystems"]] <- module@scaledOutputs[["NO2Emit"]]*getCF(data.dir,"nitrogen oxides (as nitrogen dioxide)","OzoneFormationEcosystem",perspective,CFs)
  module@foreground[["Terrestrial Acidification"]] <- module@scaledOutputs[["NO2Emit"]]*getCF(data.dir,"NOx","TerrestrialAcidification",perspective,CFs)
  module@foreground[["Fine Particulate Matter Formation"]] <- module@scaledOutputs[["NO2Emit"]]*getCF(data.dir,"NOx","ParticulateMatter",perspective,CFs)
  
  #Return object
  return(module)
}

#Define fermenter object
fermenter <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Fermenter",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["AirInput"]] <- flowDists[["BioreactorAir_(N2/O2)Input"]]
  module@parameters[["AirOutput"]] <- flowDists[["BioreactorAirOutput"]]
  module@parameters[["CO2Output"]] <- flowDists[["BioreactorCO2Output"]]
  module@parameters[["Fermenters"]] <- flowDists[["GeneralNumber_of_fermentersParameter"]]
  module@parameters[["FermPlant"]] <- flowDists[["GeneralNumber_of_fermentersParameter"]]/parameters[["BaseFermSize"]][run]/parameters[["times"]][run]
  module@parameters[["NitrogenReq"]] <- flowDists[["BioreactorNitrogenInput"]]
  module@parameters[["NitrogenInput"]] <- module@parameters[["NitrogenReq"]]*parameters[["AmmoniumNitrogen"]][run]
  module@parameters[["SaltInput"]] <- flowDists[["BioreactorSaltsInput"]]
  module@parameters[["SugarInput"]] <- flowDists[["BioreactorFeedstockInput"]]
  module@parameters[["Product_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["RecycledWater"]] <- (flowDists[["CakeWater"]]-flowDists[["EvaporatorCakeOutput"]]*parameters[["WaterContentCake"]][run])+(flowDists[["DryingWater_vaporOutput"]]-parameters[["WaterContentBiomass"]][run]*flowDists[["DryingBiomass_Output"]])
  module@parameters[["WaterReq"]] <- flowDists[["BioreactorWaterInput"]]
  module@parameters[["RunTime"]] <- flowDists[["GeneralFermentation_Run_TimeParameter"]]
  module@parameters[["FermVolume"]] <- module@parameters[["Fermenters"]] * parameters[["FermenterMass"]][run]
  module@parameters[["RunningTime"]] <- parameters[["AnnualTime"]][run]/(module@parameters[["RunTime"]]+parameters[["DownTime"]][run])*module@parameters[["RunTime"]]
  module@parameters[["WaterInput"]] <- module@parameters[["WaterReq"]] - module@parameters[["RecycledWater"]]
  #Lignocellulosic sugar arrives in solution so water input needs correcting
  if(parameters$geogs[run] == "USLigno") {
    module@parameters[["WaterInput"]] <- module@parameters[["WaterReq"]] - module@parameters[["SugarInput"]]*parameters[["NRELSugarDilution"]][run]
  }
  module@parameters[["ElecRequired"]] <- parameters[["AgitationAeration"]][run]/1000*module@parameters[["FermVolume"]]*module@parameters[["RunningTime"]]
  module@parameters[["SteamReq"]] <- parameters[["KgtoMJSteam"]][run]*parameters[["Sterilisation"]][run]*(module@parameters[["NitrogenInput"]]+module@parameters[["SaltInput"]]+
                                                                                           module@parameters[["SugarInput"]]+module@parameters[["WaterInput"]])
  module@parameters[["EmbodiedCarbon"]] <- parameters[[paste0("EmbodiedCarbon",mol)]][run]*flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["FermBroth"]] <- (flowDists[["BioreactorFeedstockInput"]] + flowDists[["BioreactorNitrogenInput"]] +
                                        flowDists[["BioreactorSaltsInput"]] + flowDists[["BioreactorWaterInput"]] +
                                        flowDists[["BioreactorAir_(N2/O2)Input"]] + flowDists[["FertiliserMiscMassInput"]] -
                                        flowDists[["BioreactorCO2Output"]] - flowDists[["BioreactorAirOutput"]])

  #Define inputs
  module@inputs[["Electricity"]] <- module@parameters[["ElecRequired"]]
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["FermPlant"]] <- module@parameters[["FermPlant"]]
  module@inputs[["Nitrogen"]] <- module@parameters[["NitrogenInput"]]
  module@inputs[["Salt"]] <- module@parameters[["SaltInput"]]
  module@inputs[["Sugar"]] <- module@parameters[["SugarInput"]]
  module@inputs[["Water"]] <- module@parameters[["WaterInput"]]
  module@inputs[["EmbodiedCarbon"]] <- module@parameters[["EmbodiedCarbon"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_out"]]
  module@outputs[["CO2Emit"]] <- module@parameters[["CO2Output"]]
  module@outputs[["Fermentation Broth"]] <- module@parameters[["FermBroth"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Electricity"]] <- scaleList(impactList(paste0("Electricity",parameters$energySource[run]),background),module@scaledInputs[["Electricity"]])
  module@background[["Nitrogen"]] <- scaleList(impactList(parameters$NitrogenSource[run],background),module@scaledInputs[["Nitrogen"]])
  #If US ligno scenario sugar impacts are accounted for elsewhere
  if(!parameters$geogs[run] == "USLigno"){
    module@background[["Sugar"]] <- scaleList(impactList("Sugar",background),module@scaledInputs[["Sugar"]])
  }
  module@background[["Water"]] <- scaleList(impactList("Water",background),module@scaledInputs[["Water"]])
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])
  module@background[["FermPlant"]] <- scaleList(impactList("EthanolFermPlant",background),module@scaledInputs[["FermPlant"]])
  module@background[["SodiumChloride"]] <- scaleList(impactList("SodiumChloride",background),module@scaledInputs[["Salt"]])
  module@background[["EmbodiedCarbon"]] <- impactList()
  #Add credit for embodied carbon
  module@background$EmbodiedCarbon$`Climate change, default, excl biogenic carbon` <- -module@scaledInputs[["EmbodiedCarbon"]]
  
  #Define foreground impacts
  module@foreground[["Climate change, incl biogenic carbon"]] <- module@scaledOutputs[["CO2Emit"]]*getCF(data.dir,"Carbon dioxide","GlobalWarming",perspective,CFs)
  
  #Return object
  return(module)
}

#Define centrifuge object
centrifuge <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Centrifuge",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["Biomass"]] <- flowDists[["DryingBiomass_Output"]]
  module@parameters[["Product_in_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["WaterAdded"]] <- flowDists[["CentrifugeWaterInput"]]
  module@parameters[["FermBroth"]] <- (flowDists[["BioreactorFeedstockInput"]] + flowDists[["BioreactorNitrogenInput"]] +
                                         flowDists[["BioreactorSaltsInput"]] + flowDists[["BioreactorWaterInput"]] +
                                         flowDists[["BioreactorAir_(N2/O2)Input"]] + flowDists[["FertiliserMiscMassInput"]] -
                                         flowDists[["BioreactorCO2Output"]] - flowDists[["BioreactorAirOutput"]])
  module@parameters[["ElecRequired"]] <- parameters[["CentrifugeEcoli"]][run]/1000*(module@parameters[["FermBroth"]]+module@parameters[["WaterAdded"]])
  module@parameters[["FermBrothIn"]] <- module@parameters[["FermBroth"]] + flowDists[["CentrifugeWaterInput"]]
  module@parameters[["FermBrothOut"]] <- module@parameters[["FermBroth"]] - (flowDists[["DryingWater_vaporOutput"]] + flowDists[["DryingBiomass_Output"]])
  module@parameters[["WetBiomass"]] <- flowDists[["DryingWater_vaporOutput"]] + flowDists[["DryingBiomass_Output"]]
  
  #Define inputs
  module@inputs[["Electricity"]] <- module@parameters[["ElecRequired"]]
  module@inputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@inputs[["Water"]] <- module@parameters[["WaterAdded"]]
  module@inputs[["Fermentation Broth"]] <- module@parameters[["FermBrothIn"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@outputs[["Biomass"]] <- module@parameters[["Biomass"]]
  module@outputs[["Fermentation Broth"]] <- module@parameters[["FermBrothOut"]]
  module@outputs[["Wet Biomass"]] <- module@parameters[["WetBiomass"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Electricity"]] <- scaleList(impactList(paste0("Electricity",parameters$energySource[run]),background),module@scaledInputs[["Electricity"]])
  module@background[["Water"]] <- scaleList(impactList("Water",background),module@scaledInputs[["Water"]])
  
  #Return object
  return(module)
}

#pH Adjustment module
phAdjustment <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="pH Adjustment",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["NaOH"]] <- flowDists[["pH_adjustmentNaOHInput"]]
  module@parameters[["Product_in_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["FermBrothIn"]] <- (flowDists[["BioreactorFeedstockInput"]] + flowDists[["BioreactorNitrogenInput"]] +
                                         flowDists[["BioreactorSaltsInput"]] + flowDists[["BioreactorWaterInput"]] +
                                         flowDists[["BioreactorAir_(N2/O2)Input"]] + flowDists[["FertiliserMiscMassInput"]] -
                                         flowDists[["BioreactorCO2Output"]] + flowDists[["CentrifugeWaterInput"]] - 
                                       (flowDists[["BioreactorAirOutput"]] + flowDists[["DryingWater_vaporOutput"]] +
                                          flowDists[["DryingBiomass_Output"]]))
  module@parameters[["FermBrothOut"]] <- module@parameters[["FermBrothIn"]] + flowDists[["pH_adjustmentNaOHInput"]]
  

  #Define inputs
  module@inputs[["NaOH"]] <- module@parameters[["NaOH"]]
  module@inputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@inputs[["Fermentation Broth"]] <- module@parameters[["FermBrothIn"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@outputs[["Fermentation Broth"]] <- module@parameters[["FermBrothOut"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["NaOH"]] <- scaleList(impactList("NeutralisingAgent",background),module@scaledInputs[["NaOH"]])
  
  #Return object
  return(module)
}

#Extracter module
extractor <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Extractor",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["ButanolReq"]] <- flowDists[["ExtracterSolventInput"]]
  module@parameters[["ButanolShift"]] <- flowDists[["ExtracterSolventInput"]]/parameters[["SolventLossRate"]][run]
  module@parameters[["ButanolRec"]] <- module@parameters[["ButanolShift"]]*(1-parameters[["SolventLossRate"]][run])
  module@parameters[["Product_in_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["CakeOutput"]] <- flowDists[["EvaporatorCakeOutput"]]
  module@parameters[["FermBrothIn"]] <- (flowDists[["EvaporatorCakeOutput"]] + flowDists[["Distillation_step_2Heavy_OrganicsOutput"]] +
                                      flowDists[["Distillation_step_2ProductOutput"]] + flowDists[["EvaporatorWaterOutput"]] +
                                      flowDists[["EvaporatorWater_vaporOutput"]]) - flowDists[["ExtracterSolventInput"]]
  module@parameters[["ExtractedProduct"]] <- (module@parameters[["ButanolRec"]] + flowDists[["Distillation_step_2Heavy_OrganicsOutput"]] +
                                       flowDists[["Distillation_step_2ProductOutput"]])
  module@parameters[["WetBiomassCake"]] <- (flowDists[["EvaporatorCakeOutput"]] + flowDists[["EvaporatorWaterOutput"]] +
                                        flowDists[["EvaporatorWater_vaporOutput"]])
                              
  #Define inputs
  module@inputs[["Butanol"]] <- module@parameters[["ButanolReq"]]
  module@inputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@inputs[["ButanolRec"]] <- module@parameters[["ButanolRec"]]
  module@inputs[["Fermentation Broth"]] <- module@parameters[["FermBrothIn"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@outputs[["Cake"]] <- module@parameters[["CakeOutput"]]
  module@outputs[["Extracted Product"]] <- module@parameters[["ExtractedProduct"]]
  module@outputs[["Wet Biomass/Cake"]] <- module@parameters[["AuxFlowOut"]]
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Butanol"]] <- scaleList(impactList("1-Butanol",background),module@scaledInputs[["Butanol"]])
  
  #Return object
  return(module)
}

#Distillation 1 module
distillation1 <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Stripper",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["Product_in_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["ButanolShift"]] <- flowDists[["ExtracterSolventInput"]]/parameters[["SolventLossRate"]][run]
  module@parameters[["ButanolRec"]] <- module@parameters[["ButanolShift"]]*(1-parameters[["SolventLossRate"]][run])
  module@parameters[["SteamReq"]] <- (module@parameters[["ButanolShift"]]*(parameters[["ButanolSensHeat"]][run]*parameters[["Dist1DeltaT"]][run] + parameters[["ButanolVapHeat"]][run]) +
                                        flowDists[["Distillation_step_2ProductOutput"]]*parameters[[paste0(mol,"SensHeat")]][run]*parameters[["Dist1DeltaT"]][run])/parameters[["DistEff"]][run]
  module@parameters[["ExtractedProduct"]] <- (module@parameters[["ButanolRec"]] + flowDists[["Distillation_step_2Heavy_OrganicsOutput"]] +
                                       flowDists[["Distillation_step_2ProductOutput"]])
  module@parameters[["ImpureProduct"]] <- (flowDists[["Distillation_step_2Heavy_OrganicsOutput"]] +
                                        flowDists[["Distillation_step_2ProductOutput"]])
  module@parameters[["Recycled Butanol"]] <- (module@parameters[["ButanolRec"]])
  
  #Define inputs
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@inputs[["Extracted Product"]] <- module@parameters[["ExtractedProduct"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@outputs[["Impure Product"]] <- module@parameters[["ImpureProduct"]]
  module@outputs[["Recycled Butanol"]] <- module@parameters[["RecycledButanol"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])

  #Return object
  return(module)
}

#Distillation 2 module
distillation2 <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Distillation",scaling="Product",foreground=impactList())
  #Define parameters
  module@parameters[["Product_in_out"]] <- flowDists[["Distillation_step_2ProductOutput"]]
  module@parameters[["SteamReq"]] <- (flowDists[["Distillation_step_2ProductOutput"]]*
                                        (parameters[[paste0(mol,"SensHeat")]][run]*
                                           parameters[[paste0("Dist2DeltaT",mol)]][run] + 
                                           parameters[[paste0(mol,"VapHeat")]][run]))/parameters[["DistEff"]][run]
  module@parameters[["Heavy_Organics"]] <- flowDists[["Distillation_step_2Heavy_OrganicsOutput"]]
  module@parameters[["ImpureProduct"]] <- flowDists[["Distillation_step_2Heavy_OrganicsOutput"]]+flowDists[["Distillation_step_2ProductOutput"]]
  
  #Define inputs
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@inputs[["Impure Product"]] <- module@parameters[["ImpureProduct"]]
  
  #Define outputs
  module@outputs[["Product"]] <- module@parameters[["Product_in_out"]]
  module@outputs[["Biomass"]] <- module@parameters[["Heavy_Organics"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Product"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])

  #Return object
  return(module)
}

#Evaporator Biomass module
evaporatorBiomass <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Evaporator (Biomass)",scaling="Biomass",foreground=impactList())
  #Define parameters
  module@parameters[["Biomass"]] <- flowDists[["DryingBiomass_Output"]]
  module@parameters[["WaterVapour"]] <- flowDists[["DryingWater_vaporOutput"]]-parameters[["WaterContentBiomass"]][run]*flowDists[["DryingBiomass_Output"]]
  module@parameters[["SteamReq"]] <- module@parameters[["WaterVapour"]]*parameters[["EvaporationTripleEffect"]][run]*parameters[["KgtoMJSteam"]][run]
  module@parameters[["WetBiomassIn"]] <- module@parameters[["Biomass"]] + module@parameters[["WaterVapour"]]
  module@parameters[["BiomassOut"]] <- module@parameters[["Biomass"]] + parameters[["WaterContentBiomass"]][run]*flowDists[["DryingBiomass_Output"]]
  module@parameters[["WaterEvap"]] <- module@parameters[["WaterVapour"]] 
  
  #Define inputs
  module@inputs[["Biomass"]] <- module@parameters[["Biomass"]]
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["Wet Biomass"]] <- module@parameters[["WetBiomassIn"]]
  
  #Define outputs
  module@outputs[["Biomass"]] <- module@parameters[["Biomass"]]
  module@outputs[["WaterVapour"]] <- module@parameters[["WaterVapour"]]
  module@outputs[["Biomass"]] <- module@parameters[["BiomassOut"]]
  module@outputs[["Recycled Water"]] <- module@parameters[["WaterEvap"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Biomass"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])

  #Return object
  return(module)
}

#Evapourator module - Cake to dryer and waste
evaporatorCake <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Evaporator (Cake)",scaling="Cake",foreground=impactList())
  #Define parameters
  module@parameters[["WaterVapour"]] <- flowDists[["EvaporatorCakeOutput"]]*parameters[["WaterContentCake"]][run]
  module@parameters[["WaterEvap"]] <- flowDists[["CakeWater"]] - module@parameters[["WaterVapour"]]
  module@parameters[["Cake"]] <- flowDists[["EvaporatorCakeOutput"]]
  module@parameters[["SteamReq"]] <- module@parameters[["WaterEvap"]]*parameters[["EvaporationTripleEffect"]][run]*parameters[["KgtoMJSteam"]][run]
  module@parameters[["WetCakeIn"]] <- (flowDists[["EvaporatorCakeOutput"]] + flowDists[["CakeWater"]])
  module@parameters[["WetCakeOut"]] <- flowDists[["EvaporatorCakeOutput"]] + module@parameters[["WaterVapour"]]
  module@parameters[["WaterEvap"]] <- module@parameters[["WaterEvap"]]
  
  #Define inputs
  module@inputs[["Cake"]] <- module@parameters[["Cake"]]
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["Wet Cake/Biomass"]] <- module@parameters[["WetCake"]]
  
  #Define outputs
  module@outputs[["Cake"]] <- module@parameters[["Cake"]]
  module@outputs[["Wet Cake/Biomass"]] <- module@parameters[["WetCakeOut"]]
  module@outputs[["Recycled Water"]] <- module@parameters[["WaterEvap"]]
  
  #Scale results
  scalingFactor <- requirement/module@outputs[["Cake"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])
  
  #Return object
  return(module)
}

#Dryer module
dryer <- function(data.dir,requirement = 1, run, parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Dryer",scaling="Cake",foreground=impactList())
  #Define parameters
  module@parameters[["WaterVapour"]] <- flowDists[["EvaporatorCakeOutput"]]*parameters[["WaterContentCake"]][run]
  module@parameters[["WaterEvap"]] <- module@parameters[["WaterVapour"]] - flowDists[["EvaporatorCakeOutput"]]*parameters[["WaterContentWaste"]][run]
  module@parameters[["Cake"]] <- flowDists[["EvaporatorCakeOutput"]]
  module@parameters[["SteamReq"]] <- module@parameters[["WaterEvap"]]*parameters[["DryingSteam"]][run]*parameters[["KgtoMJSteam"]][run]
  module@parameters[["WetCake"]] <- module@parameters[["WaterVapour"]] + flowDists[["EvaporatorCakeOutput"]]
  
  #Define inputs
  module@inputs[["Cake"]] <- module@parameters[["Cake"]]
  module@inputs[["Steam"]] <- module@parameters[["SteamReq"]]
  module@inputs[["Wet Cake"]] <- module@parameters[["WetCake"]]
  
  #Define outputs
  module@outputs[["Waste"]] <- module@parameters[["Cake"]]
  module@outputs[["WaterVapour"]] <- module@parameters[["WaterVapour"]]
  
  #Scale results
  scalingFactor <- requirement/module@inputs[["Cake"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),module@scaledInputs[["Steam"]])
  module@background[["Waste"]] <- scaleList(impactList("WasteTreatment2",background),module@scaledOutputs[["Waste"]])

  #Return object
  return(module)
}

#Combustor module
combustor <- function(data.dir,requirement = 1,run,parameters,flowDists,background,CFs,perspective,mol="Cadaverine") {
  #Define parameters
  module <- new("ProcessModule",name="Combustor",scaling="Biomass",foreground=impactList())
  #Define parameters
  module@parameters[["Biomass"]] <- 1
  module@parameters[["SteamGen"]] <- module@parameters[["Biomass"]]*parameters[["BiomasstoHeat"]][run]
  module@parameters[["CombustorUsage"]] <- module@parameters[["Biomass"]]/parameters[["BiomassPerMJCombustor"]][run]
  #Define inputs
  module@inputs[["Biomass"]] <- module@parameters[["Biomass"]]
  module@inputs[["CombustorUsage"]] <- module@parameters[["CombustorUsage"]]
  
  #Define outputs
  module@outputs[["Steam"]] <- module@parameters[["SteamGen"]]
  
  #Scale results
  scalingFactor <- requirement/module@inputs[["Biomass"]]
  module@scaledOutputs <- scaleList(module@outputs,scalingFactor)
  module@scaledInputs <- scaleList(module@inputs,scalingFactor)
  
  #Define background impacts
  module@background[["Steam"]] <- scaleList(impactList(paste0("Steam",parameters$energySource[run]),background),-module@scaledOutputs[["Steam"]])
  module@background[["Combustor"]] <- scaleList(impactList("Combustor",background),module@scaledInputs[["CombustorUsage"]])

  #Return object
  return(module)
}