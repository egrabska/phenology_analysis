#modelling LSP from SENTINEL-2-----------------------------------

packages_list = c("raster", "sf", "exactextractr", "dplyr", "ggplot2", "tidyverse", "tsibble", "fable", "bfast",
                  "matrixStats", "signal", "data.table", "ggpubr", "fabletools",
                  "plotly", "mgcv", "phenopix", "MASS", "expss", "data.table", "forecast")
lapply(packages_list, require, character.only = TRUE)

#1--------Cleaning tables in .csv format acquired from GEE-------------

#df in csv format containing variables in this order: date, index, sample ID, species
df =  read.csv("C:/Ewa/Barlinek_reviews2/Barlinek_2017_2022_MTCI_pts.csv")

df$system.index = substr(df$system.index, 1, 8) %>%
  as.Date(format =  "%Y%m%d")
names(df) = c("date", "index", "id", "species")

#filtering - removing NA and calculating means for duplicated values (e.g. from two different S2 tiles)
df2 = df %>% 
  drop_na() %>%
  group_by(date, species, id) %>%
  summarise(index = mean(index)) %>% 
  ungroup() 

unique(df2$date)

summary(df2)
#removing low summer values
df_sum = df2 %>%
 dplyr::filter(date >= as.Date(paste(year(date), 05, 15, sep = "-")),
               date <= as.Date(paste(year(date), 09, 01, sep = "-")))

df_win = df2 %>%
 dplyr::filter(date < as.Date(paste(year(date), 05, 15, sep = "-")) |
          date > as.Date(paste(year(date), 09, 01, sep = "-")))

df_sum[, 4][df_sum[, 4] < 2] = NA #mtci 2 evi 0.3 ndvi 0.35
df2 = rbind(df_sum, df_win)

#removing highest and lowest values and December-February observations
df3 = df2 %>% 
  dplyr::filter(index > 0 & index < 8) %>% #in case of MTCI > 0 & MTCI < 8 NDMI -0.2;0.65
  dplyr::filter(!grepl("-01-", date)) %>%
  dplyr::filter(!grepl("-02-", date)) %>%
  arrange(date)


 
#2----modelling time series and detecting SOS/EOS based on derivatives
#select single year observations
sel = df3 %>%
  dplyr::filter(date > "2018-03-10" & 
           date < "2018-11-20")  

#check visualization
ggplot(sel, aes(date, index, colour = species))+
  geom_point(size = 2, alpha = 0.4)


#function to model time series using GAM and detecting SOS/EOS dates based on derivatives
gam_multiple_GEE = function(id_no, input_df){
  options(warn=-1)
  df_long = input_df %>%
    dplyr::filter(id == id_no) 
  df_ts =df_long[,4] %>%
    ts() %>%
    tsclean(iterate = 2) %>%
    as.data.frame()
  df_ts$date = df_long$date
  df_bfast = bfastts(df_ts$index, df_ts$date, type = c("irregular", "1-day"))
  df_tibble = tibble(date = seq(as.Date(df_long$date[1]), by = "day", 
               length.out = length(df_bfast)), value = df_bfast) %>%
    as_tsibble(index = date) %>%
    fill_gaps() %>%
    ts() %>% 
    as.data.frame()
  model = gamm(df_tibble[,2] ~ s(date, bs = "cr", k = 12), 
               data = df_tibble, method = "REML")
  df_tibble$predicted = predict.gam(model$gam, df_tibble)
  df_tibble$date = anytime::anydate(df_tibble$date)
  df_ts2 = ts(df_tibble) %>%
    as.data.frame()
  derivative = diff(df_ts2$predicted)/diff(df_ts2$date)
  der = cbind(df_tibble[-3,], round(derivative,4))
  der = der[,c(-2,-3)]
  names(der) = c("date", paste(id_no))
  der_long = der %>%
    gather("id", "value", 2:length(der))
  sos = der_long %>% 
    slice(which.max(value))
  eos = der_long %>% 
    slice(which.min(value))
  der_dates = left_join(sos, eos, "id") %>%
    dplyr::select(date.x, id, date.y)
  names(der_dates) = c("SOS", "id", "EOS")
  der_dates[,c('SOS', 'EOS')] = sapply(der_dates[,c('SOS', 'EOS')],strftime, format = "%j")
  return(der_dates)
}

#how it works for single pixel, e.g. id =1
single = gam_multiple_GEE(1, sel)
single

#create list of unique ids(for single id and single date more than 15 observations)
obs_n = sel %>%
  dplyr::count(id) %>%
  dplyr::filter(n > 15)

lista = unique(obs_n$id) 

#use lapply to calculate lsp for all pixels id
start = Sys.time()
multiple = lapply(lista, gam_multiple_GEE, sel) %>%
  rbindlist()
end = Sys.time()
end-start

codes = sel %>%
  dplyr::select(species, id) %>%
  distinct()
codes$id = as.character(codes$id)

#for time series analysis
multiple = left_join(mulitple, codes, "id")


ggplot(multiple, aes(as.numeric(SOS), as.numeric(EOS), color = species))+
  geom_point()+
  xlim(100, 140)+
  ylim(200, 300)


write.csv2(multiple, "C:/Ewa/Barlinek_reviews2/EVI_2022_SOS_EOS.csv")






