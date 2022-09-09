#modelling LSP from SENTINEL-2-----------------------------------

packages_list = c("raster", "sf", "exactextractr", "dplyr", "ggplot2", "tidyverse", "tsibble", "fable", "bfast",
                  "matrixStats", "signal", "data.table", "ggpubr", "fabletools",
                  "plotly", "mgcv", "phenopix", "MASS", "expss", "data.table", "forecast")
lapply(packages_list, require, character.only = TRUE)

#1--------Cleaning tables in .csv format acquired from GEE-------------

#df in csv format containign variables (in this order): date, index, sample ID, species
df =  read.csv("C:/Ewa/Barlinek_reviews2/Barlinek_2017_2022_EVI.csv") 

df$system.index = substr(df$system.index, 1, 8) %>%
  as.Date(format =  "%Y%m%d")
names(df) = c("date", "index", "id", "species")

#filtering - removing NA and calculating means for duplicated values (e.g. from two different S2 tiles)
df2 = df %>% 
  drop_na() %>%
  group_by(date, species, id) %>%
  dplyr::summarise_all(funs(mean)) %>% 
  ungroup() 


#removing highest and lowest values and December-February observations
df3 = df2 %>% 
  filter(index < 1.2 & index> 0) %>% #in case of MTCI > 0 & MTCI < 8
  filter(!grepl("-01-", date)) %>%
  filter(!grepl("-12-", date)) %>%
  filter(!grepl("-02-", date)) %>%
  arrange(date)

#check visualization
ggplot(df3, aes(date, index, colour = species))+
  geom_point(size = 2, alpha = 0.4)


#2----modelling time series and detecting SOS/EOS based on derivatives
#select single year observations
sel = df3 %>%
  filter(date > "2019-03-15" & 
           date < "2019-11-30")  



#function to model time series using GAM and detecting SOS/EOS dates based on derivatives
gam_multiple_GEE = function(id_no, input_df){
  options(warn=-1)
  df_long = input_df %>%
    dplyr::filter(id == id_no) 
  df_ts =df_long[,4] %>%
    ts() %>%
    tsclean() %>%
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
  #only smoothing/modelling------------------------
  #sp.tibble = sp.tibble[,c(1,3)]
  #names(sp.tibble) = c("date", paste(id))
  #sp.tibble_long = sp.tibble %>%
  #gather("id", "value", 2:length(sp.tibble))
  #return(sp.tibble_long)
  #---end------------------------------------------
  df_ts2 = ts(df_tibble) %>%
    as.data.frame()
  derivative = diff(df_ts2$predicted)/diff(df_ts2$date)
  der = cbind(df_tibble[-3,], round(derivative,4))
  der = der[,c(-2,-3)]
  names(der) = c("date", paste(id))
  der_long = der %>%
    gather("id", "value", 2:length(der))
  sos = der_long %>% 
    slice(which.max(value))
  eos = der_long %>% 
    slice(which.min(value))
  der_dates = left_join(sos, eos, "id") %>%
    dplyr::select(date.x, id, date.y)
  names(der_dates) = c("SOS", "id", "EOS")
  der_dates$SOS = strftime(der_dates$SOS, format = "%j")
  der_dates$EOS = strftime(der_dates$EOS, format = "%j")
  return(der_dates)
}


single = gam_multiple_GEE(10, sel)


obs_n = sel %>%
  dplyr::count(id) %>%
  dplyr::filter(n > 20)

#create list of unique ids(for single id and single date more than 20 observations)
lista = unique(obs_n$id) 

start = Sys.time()
mulitple = lapply(lista[1:50], gam_multiple_GEE, sel) %>%
  rbindlist()
end = Sys.time()
end-start

species_id = sel[,2:3]
kody = species_id %>%
  distinct()
kody$id = as.character(kody$id)

#for time series analysis
mulitple = left_join(mulitple, kody, "id")
str(mulitple)

ggplot(mulitple, aes(as.numeric(SOS), as.numeric(EOS), color = species))+
  geom_point()+
  xlim(110, 150)+
  ylim(270, 300)


write.csv2(sel_gam, "C:/Ewa/Barlinek_reviews2/MTCI_2019_SOS_EOS.csv")


id_no = 10
input_df =sel


ggplot(sp_ind_long, aes(name, value))+
  geom_point()
names(sp_ind_long) = c("date", "species", "id", "value_old")
sp_ind_long$date = as.Date(sp_ind_long$date, format = "%Y-%m-%d")

sp_final = left_join(sp.tibble, sp_ind_long, "date")

ggplot(sp_final)+
  #geom_point(aes(date, value_old), color = "red", size = 1.9, alpha = 1)+
  #geom_point(aes(date, value), color = "black", size = 1.9, alpha = 1)+
  geom_line(aes(date, predicted), color = "darkgreen", size = 0.8)+
  geom_vline(xintercept = as.numeric(as.Date("2019-04-22")), linetype="dashed", 
             color = "black", size = 1.0)+
  geom_vline(xintercept = as.numeric(as.Date("2019-08-20")), linetype="dashed", 
             color = "black", size = 1.0)+
  theme_bw()



#applying double log + derivatives--------------------------------------------
x = tab %>%
  dplyr::count(species_cd)
x
tab_sp = tab %>% 
  filter(species_cd == "TP")

start = Sys.time()
df_total = data.frame()
for (i in 1:369){
  model = sg_multiple_GEE(i, tab_sp) #version with spread
  df = data.frame(model)
  df_total = rbind(df_total,df)
}
end = Sys.time()
end - start

cleaned_sp$id = seq.int(nrow(cleaned_sp)) %>%
  as.character()
cleaned_date = left_join(df_total, cleaned_sp[,c('a_i_num', 'id')], "id")
cleaned_date$SOS = strftime(cleaned_date$SOS, format = "%j") %>%
  as.numeric()
cleaned_date$EOS = strftime(cleaned_date$EOS, format = "%j") %>%
  as.numeric()

cleaned_date = cleaned_date %>% 
  filter(EOS > 150 & EOS < 300 & SOS < 200)

ggplot(cleaned_date, aes(SOS, EOS))+
  geom_point()


#joining with environmental data--------------------------------
envi = read.csv2("D:/18_Phenology/zachpom/envi.csv")
envi$a_i_num = as.character(envi$a_i_num)
cleaned_date$a_i_num = as.character(cleaned_date$a_i_num)
dates = left_join(cleaned_date, envi, "a_i_num")
str(dates)

write.csv2(dates, "D:/18_Phenology/GEE/DB_MTCI2021.csv")

str(dates)
p = ggplot(dates, aes(age, SOS, color = water_distance.tif))+
  geom_point()
#ylim(250, 290)

ggplotly(p)

m1 = gam(SOS ~ s(age) + s(eudem_dem_5deg_combinedPL_92_exp.tif) + s(built_up50.tif) +
           s(water_distance.tif), data = dates)
summary(m1)


input_df = tab

dl_multiple_GEE = function(id, input_df) {
  sp_ind = input_df %>%
    dplyr::select(starts_with("2"))
  sp_ind$ID = seq.int(nrow(sp_ind))
  for (i in 1:(length(sp_ind) - 1)){              
    x = sp_ind[i]
    q1 = quantile(x, na.rm = TRUE)[2]
    q3 = quantile(x, na.rm = TRUE)[4]
    iqr = q3-q1
    x[x < q1 - iqr*1.5 | x > q3 + iqr*1.5] = NA
    sp_ind[,i] = x
  }
  sp_ind_long = tidyr::pivot_longer(sp_ind, cols = 1:length(sp_ind) - 1) %>%
    group_by(ID)
  sp_ind_long$name = as.Date(sp_ind_long$name, format =  "%Y-%m-%d") 
  sel = sp_ind_long %>%
    dplyr::filter(ID== id)
  sp.bfast = bfastts(sel$value, sel$name, type = c("irregular", "1-day"))
  sp.tibble = tibble(
    date = seq(as.Date(sel$name[1]), by = "day", length.out = length(sp.bfast)), 
    value = sp.bfast
  ) %>%
    as_tsibble(index = date) %>%
    fill_gaps() 
  sp.inter = sp.tibble %>%
    model(naive = ARIMA(value ~ -1 + pdq(0,1,0) + PDQ(0,0,0))) %>%
    interpolate(sp.tibble)
  double_log2018 = sp.inter[sp.inter$date >= "2021-04-01" & sp.inter$date <= "2021-10-11",] %>%
    as.ts(frequency = 365) %>%
    FitDoubleLogBeck()
  pred2018 = double_log2018$predicted %>%
    as.data.frame()
  pred2018$date = seq(as.Date("2018/04/01"), as.Date("2018/10/31"), "day")
  names(pred2018) = c(paste(id), "date")
  y.ts = ts(pred2018) %>%
    as.data.frame()
  derivative = diff(y.ts[,1])/diff(y.ts$date)
  der = cbind(pred2018[-1,], round(derivative,4))
  der = der[,-1]
  names(der) = c("date", paste(id))
  der_long = der %>%
    gather("id", "value", 2:length(der))
  sos = der_long %>% 
    slice(which.max(value))
  eos = der_long %>% 
    slice(which.min(value))
  der_dates = left_join(sos, eos, "id") %>%
    dplyr::select(date.x, id, date.y)
  names(der_dates) = c("SOS", "id", "EOS")
  return(der_dates)
}


sp_ind$ID

id = 18 
input_df = tab
lista = unique(input_df$a_i_num)
lista[1]
id = 160801119
#what about sg_multiple_GEE
sg_multiple_GEE = function(id, input_df){
  sp_ind = input_df[,3:length(input_df)]
  for (i in 2:(length(sp_ind))){              
    x = sp_ind[i]
    q1 = quantile(x, na.rm = TRUE)[2]
    q3 = quantile(x, na.rm = TRUE)[4]
    iqr = q3-q1
    x[x < q1 - iqr*1.5 | x > q3 + iqr*1.5] = NA
    sp_ind[,i] = x
  }
  sp_ind_long = tidyr::pivot_longer(sp_ind, cols = 2:length(sp_ind)) %>%
    group_by(a_i_num)
  sp_ind_long$name = as.Date(sp_ind_long$name, format =  "%Y-%m-%d") 
  sel = sp_ind_long %>%
    dplyr::filter(a_i_num== id)%>%
    arrange(name) %>% 
    drop_na()
  sp.bfast = bfastts(sel$value, sel$name, type = c("irregular", "1-day"))
  sp.tibble = tibble(
    date = seq(as.Date(sel$name[1]), by = "day", length.out = length(sp.bfast)), 
    value = sp.bfast
  ) %>%
    as_tsibble(index = date) %>%
    fill_gaps() 
  sp.inter = sp.tibble %>%
    model(naive = ARIMA(value ~ -1 + pdq(0,1,0) + PDQ(0,0,0))) %>%
    interpolate(sp.tibble)
  sg = sgolayfilt(sp.inter$value, p = 2, n = 75) %>% as.data.frame()
  sg$date = sp.inter$date
  names(sg) = c(paste(id), "date")
  y.ts = ts(sg) %>%
    as.data.frame()
  derivative = diff(y.ts[,1])/diff(y.ts$date)
  der = cbind(sg[-1,], round(derivative,4))
  der = der[,-1]
  names(der) = c("date", paste(id))
  der_long = der %>%
    gather("id", "value", 2:length(der))
  sos = der_long %>% 
    slice(which.max(value))
  eos = der_long %>% 
    slice(which.min(value))
  der_dates = left_join(sos, eos, "id") %>%
    dplyr::select(date.x, id, date.y)
  names(der_dates) = c("SOS", "id", "EOS")
  return(der_dates)
}

x = gam_multiple_GEE(160801119, tab)

gam_multiple_GEE = function(id, input_df){
  options(warn=-1)
  sp_ind = input_df[,3:length(input_df)]
  for (i in 2:(length(sp_ind))){              
    x = sp_ind[i]
    q1 = quantile(x, na.rm = TRUE)[2]
    q3 = quantile(x, na.rm = TRUE)[4]
    iqr = q3-q1
    x[x < q1 - iqr*1.5 | x > q3 + iqr*1.5] = NA
    sp_ind[,i] = x
  }
  sp_ind_long = tidyr::pivot_longer(sp_ind, cols = 2:length(sp_ind)) %>%
    group_by(a_i_num)
  sp_ind_long$name = as.Date(sp_ind_long$name, format =  "%Y-%m-%d") 
  sel = sp_ind_long %>%
    dplyr::filter(a_i_num== id)%>%
    arrange(name) %>% 
    drop_na()
  sp.bfast = bfastts(sel$value, sel$name, type = c("irregular", "1-day"))
  sp.tibble = tibble(
    date = seq(as.Date(sel$name[1]), by = "day", length.out = length(sp.bfast)), 
    value = sp.bfast) %>%
    as_tsibble(index = date) %>%
    fill_gaps() %>%
    ts() %>% 
    as.data.frame()
  model = gamm(sp.tibble[,2] ~ s(date, bs = "cr", k = 12), data = sp.tibble, method = "REML")
  print(summary(model$gam))
  sp.tibble$predicted = predict.gam(model$gam, sp.tibble)
  sp.tibble$date = anytime::anydate(sp.tibble$date)
  y.ts = ts(sp.tibble) %>%
    as.data.frame()
  derivative = diff(y.ts[,3])/diff(y.ts$date)
  der = cbind(sp.tibble[-3,], round(derivative,4))
  der = der[,c(-2,-3)]
  names(der) = c("date", paste(id))
  der_long = der %>%
    gather("id", "value", 2:length(der))
  sos = der_long %>% 
    slice(which.max(value))
  eos = der_long %>% 
    slice(which.min(value))
  der_dates = left_join(sos, eos, "id") %>%
    dplyr::select(date.x, id, date.y)
  names(der_dates) = c("SOS", "id", "EOS")
  return(der_dates)
}




x = cbind(sp.inter, sg[,1])
names(sel) = c("ID", "date", "value_old")
y = left_join(x, sel, "date")
yy = left_join(y, der, "date")
names(yy) = c("date", "value",  "pred2018","ID",  "value_old", "deriv")


ggplot(yy)+
  geom_point(aes(date, value_old), color = "red")+
  geom_line(aes(date, pred2018), size = 0.9)+
  geom_line(aes(date, deriv*10 + 0.6), linetype=2)+ #on second y scale
  geom_vline(xintercept=as.numeric(der_dates$SOS), linetype=1)+
  geom_vline(xintercept=as.numeric(der_dates$EOS), linetype=1)+
  scale_y_continuous(name="index", sec.axis = sec_axis(~ . /10 - 0.06, name="derivative"))+
  theme_bw()


