# Code used to generate microclimate data and body temperature simulations of the paper 
# "Body temperature and activity patterns modulate glucocorticoid levels across lizard species: A macrophysiological approach"
# Juan G. Rubalcaba & Blanca Jimeno
# Frontiers in Ecology and Evolution, 2022
# DOI: 10.3389/fevo.2022.1032083

dir <- "..." # Insert working directory

###### Hormone data: HormoneBase ###### 

load(file=paste0(dir,"GCdata_reptiles_Hormonebase_Meiri2018_merged.RData"))

###### Extracting microclimatic data ###### 

require(NicheMapR)
require(dplyr)
require(tidyr)
require(lubridate)

i = 1 # Recommend extracting microclimate data manualy (not using "for" loops) for each observation to avoid errors
GCdata_reptiles$Binomial[i]
  
# check dates for each observation
longlat <- as.numeric(c(x=GCdata_reptiles$Longitude[i], y=GCdata_reptiles$Latitude[i]))
year1 <- as.numeric(GCdata_reptiles$Year_1[i])
year2 <- as.numeric(GCdata_reptiles$Year_final[i])
years <- year1:year2
months <- which(GCdata_reptiles[i,19:30] == "Y")
dfinish <- paste0(day,"/0",months[length(months)],"/",years[length(years)])

years
months

microclimate <- micro_ncep(loc = longlat, IUV = 0, run.gads = 0,
                           dstart = dstart, dfinish = dfinish,
                           Usrhyt = 0.01, minshade = 0, maxshade = 90)

metout <- as.data.frame(microclimate$metout)
shadmet <- as.data.frame(microclimate$shadmet)
soil <- as.data.frame(microclimate$soil)
shadsoil <- as.data.frame(microclimate$shadsoil)

microclim_data_sp <- data.frame(
  DOY = metout$DOY,
  TIME = metout$TIME,
  dates = microclimate$dates,
  SOLRAD = metout$SOLR,
  TALOC_sun = metout$TALOC,
  TSKYC_sun = metout$TSKYC,
  RHLOC_sun = metout$RHLOC,
  VLOC_sun = metout$VLOC, 
  SOILT_sun = soil$D0cm,
  TALOC_shade = shadmet$TALOC,
  TSKYC_shade = shadmet$TSKYC,
  RHLOC_shade = shadmet$RHLOC,
  VLOC_shade = shadmet$VLOC, 
  SOILT_shade = shadsoil$D0cm
)

microclim_data_sp <- microclim_data_sp %>% filter(month(microclim_data_sp$dates) %in% months)

plot(microclim_data_sp$TALOC_sun, pch=20, col="gold")
points(microclim_data_sp$TALOC_shade, pch=20, col="grey")

# Save "microclim_data_sp" as a file:

save(microclim_data_sp, file=paste0(dir,"Data/Microclimatic_data/","microclim_data_",i,".RData"))

# Merge microclim_data_sp files in a single list:

files <- list.files(paste0(dir,"Data/Microclimatic_data/"), pattern = "microclim_data")
files_numb <- as.numeric(gsub(".RData","",gsub("microclim_data_","",files)))
files_numb <- files_numb[order(files_numb)]
microclim_data <- list()
for(i in 1:length(files)){
  load(file=paste0(dir,"Data/Microclimatic_data/",files[i]))
  microclim_data[[files_numb[i]]] <- microclim_data_sp
  print(i/length(files))
}

###### Microclimate / Tb data simulations

# Function to compute operative temperature
compute_Te <- function(M, # body mass (g)
                       microclimate) # Microclimate database (requires a dataframe with the structure of "microclim_data_sp") 
  { 
  A <- 0.0314 * pi * (M/1000)^(2/3) # skin surface area (O'Connor 1999)
  l <- 3.3 * (M/1e6)^(1/3) # Body length (Mitchell 1976)
  
  # Convection heat transfer coefficient (Mitchell 1976)
  nu_sun = -1.1555e-14*(microclimate$TALOC_sun+273)^3 + 9.5728e-11*(microclimate$TALOC_sun+273)^2 + 3.7604e-08*(microclimate$TALOC_sun+273) - 3.4484e-06
  kf_sun = 1.5207e-11*(microclimate$TALOC_sun+273)^3 - 4.8574e-08*(microclimate$TALOC_sun+273)^2 + 1.0184e-04*(microclimate$TALOC_sun+273) - 3.9333e-04
  Re = l * microclimate$VLOC_sun / nu_sun  
  Nu = 0.1 * Re^0.74
  hc_sun = Nu * kf_sun / l
  
  nu_shade = -1.1555e-14*(microclimate$TALOC_shade+273)^3 + 9.5728e-11*(microclimate$TALOC_shade+273)^2 + 3.7604e-08*(microclimate$TALOC_shade+273) - 3.4484e-06
  kf_shade = 1.5207e-11*(microclimate$TALOC_shade+273)^3 - 4.8574e-08*(microclimate$TALOC_shade+273)^2 + 1.0184e-04*(microclimate$TALOC_shade+273) - 3.9333e-04
  Re = l * microclimate$VLOC_shade / nu_shade  
  Nu = 0.1 * Re^0.74
  hc_shade = Nu * kf_shade / l
  
  # Conduction heat transfer coefficient (Stevenson 1985)
  ksub = 0.027  
  ts = 0.025 * (0.001 * M / (3.1416 * 1000))^0.2 
  hg = ksub / ts
  
  # Radiative heat transfer coefficients
  Ra_sun = 4 * 0.95 * 5.670367e-8 * (microclimate$TSKYC_sun+273)^3
  Rg_sun = 4 * 0.95 * 5.670367e-8 * (microclimate$SOILT_sun+273)^3

  Ra_shade = 4 * 0.95 * 5.670367e-8 * (microclimate$TSKYC_shade+273)^3
  Rg_shade = 4 * 0.95 * 5.670367e-8 * (microclimate$SOILT_shade+273)^3
  
  # Constant j
  Ad = A * (1-0.6) 
  Ag = A * 0.6
  C = 3.7 * M
  
  a = 0.9 # absorbance to solar radiation (Gates 1980)
  j_sun = Ad / C * (a * microclimate$SOLRAD + Ra_sun * microclimate$TSKYC_sun + hc_sun * microclimate$TALOC_sun) + Ag / C * microclimate$SOILT_sun * (Rg_sun + hg)
  j_shade = Ad / C * (a * microclimate$SOLRAD + Ra_shade * microclimate$TSKYC_shade + hc_shade * microclimate$TALOC_shade) + Ag / C * microclimate$SOILT_shade * (Rg_shade + hg)
  
  # constant theta
  theta_sun <- A / C * ((1-0.6) * (Ra_sun + hc_sun) + 0.6 * (Rg_sun + hg))
  theta_shade <- A / C * ((1-0.6) * (Ra_shade + hc_shade) + 0.6 * (Rg_shade + hg))
  
  Te_sun <- j_sun / theta_sun
  Te_shade <- j_shade / theta_shade
  return(cbind(Te_shade,Te_sun))
}

Tes <- Tas <- as.data.frame(array(NA, dim=c(nrow(GCdata_reptiles),4)))
TavsTb <- TevsTb <- numeric(nrow(GCdata_reptiles))
for(i in 1:nrow(GCdata_reptiles)){
  # Air temperature
  Tas[i,1] <- mean(microclim_data[[i]]$TALOC_sun, na.rm = T) # mean Ta sun
  Tas[i,2] <- mean(microclim_data[[i]]$TALOC_shade, na.rm = T) # mean Ta shade
  Tas[i,3] <- quantile(microclim_data[[i]]$TALOC_sun, 0.95, na.rm = T) # 95%Ta sun
  Tas[i,4] <- quantile(microclim_data[[i]]$TALOC_sun, 0.05, na.rm = T) # 5% Ta sun
  
  # Operative temperature 
  M <- 10^GCdata_reptiles$log10Mass[i]
  if(!is.na(M)){
    Te <- compute_Te(M=M, microclimate = microclim_data[[i]])
    Tes[i,1] <- mean(Te[,2], na.rm = T) # mean Te sun
    Tes[i,2] <- mean(Te[,1], na.rm = T) # mean Te shade
    Tes[i,3] <- quantile(Te[,2], 0.95, na.rm = T) # 95% quart Te sun
    Tes[i,4] <- quantile(Te[,2], 0.05, na.rm = T) # 5% quart Te sun
  }
  
  # Deviation from mean species' Tb
  Tpref <- as.numeric(GCdata_reptiles$Tb)[i]
  if(!is.na(Tpref)){
    
    TavsTb[i] <- mean(abs(Tpref - rowMeans(cbind(microclim_data[[i]]$TALOC_sun, microclim_data[[i]]$TALOC_shade))))
    
    M <- 10^GCdata_reptiles$log10Mass[i]
    if(!is.na(M)){
      Te <- compute_Te(M=M, microclimate = microclim_data[[i]])
      TevsTb[i] <- mean(abs(Tpref - rowMeans(Te)))
      
    }
  }
  print(i/nrow(GCdata_reptiles))
}
