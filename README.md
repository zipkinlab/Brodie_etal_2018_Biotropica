# [Models for assessing local-scale co-abundance of animal species while accounting for differential detectability and varied responses to the environment](https://onlinelibrary.wiley.com/doi/full/10.1111/btp.12500)

### Jedediah F. Brodie, Olga E. Helmy, Jayasilan Mohd-Azlan, Alys Granados, Henry Bernard, Anthony J. Giordano, Elise Zipkin

### Biotropica

### Please contact the first author for questions about the code or data: Jedediah Brodie (jedediah.brodie@mso.umt.edu)

## Abstract
We developed a new modeling framework to assess how the local abundance of one species influences the local abundance of a potential competitor while explicitly accounting for differential responses to environmental conditions. Our models also incorporate imperfect detection as well as abundance estimation error for both species. As a case study, we applied the model to four pairs of mammal species in Borneo, surveyed by extensive and spatially widespread camera trapping. We detected different responses to elevation gradients within civet, macaque, and muntjac deer species pairs. Muntjac and porcupine species varied in their response to terrain ruggedness, and the two muntjac responded different to river proximity. Bornean endemic species of civet and muntjac were more sensitive than their widespread counterparts to habitat disturbance (selective logging). Local abundance within several species pairs was positively correlated, but this is likely due to the species having similar responses to (un-modeled) environmental conditions or resources rather than representing facilitation. After accounting for environment and correcting for false absences in detection, negative correlations in local abundance appear rare in tropical mammals. Direct competition may be weak in these species, possibly because the “ghost of competition past” or habitat filtering have already driven separation of the species in niche space. The analytical framework presented here could increase basic understanding of how ecological interactions shape patterns of abundance across the landscape for a range of taxa, and also provide a powerful tool for forecasting the impacts of global change. 

## **Data:**
There are a series of data files associated with this analysis: 1) files beginning with species names contain the species-specific count data; 2) files beginning with "sitecovs" contain covariates for the state process (those effecting species abundances); 3) files beginning with "sampcovs" contain covariates for the observation process (those effecting species detection probabilities) 


## **Code:**
[two-spp_model.R](https://github.com/zipkinlab/Brodie_etal_inpress_Biotropica/blob/master/two-spp_model.R) - R code to run the two species interaction models including code to format the data, create the BUGS text file, and run models in WinBUGS.
