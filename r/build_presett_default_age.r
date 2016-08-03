library(neotoma)
library(ggplot2)
library(maps)

source('r/utils/compile_lists_v2.r')
# source('eda/r/utils/compile_list_stepps.r')
source('r/utils/make_data_funs.r')

"%w/o%" <- function(x, y) !x %in% y #--  x without y 

# apped at the end of the filenames, especially for testing
version = 'v3'
list_name = 'kujawa'

suffix = paste0(version, '_', list_name)

###########################################################################################################
# use neotoma package to pull UMW data; save raw and then aggregate counts according to specified list_name
###########################################################################################################

# get the metadata
pollen_meta_neo  <- get_neo()
pollen_meta_neo  <- pollen_meta_neo[with(pollen_meta_neo, order(handle)),]

# write.table(pollen_meta_neo, paste('data/pollen/pollen_meta_', Sys.Date(), '.csv', sep=''), row.names=FALSE, sep=',', na='')

# download the data
pollen2k = list()
for (i in 1:length(ids)){  
  print(i)
  pollen2k[[i]] = get_download(pollen_meta_neo$datasetID[i])
  pollen2k[[i]]$dataset$state = pollen_meta_neo$state[i]
  pollen2k[[i]]$dataset$elev = pollen_meta_neo$elev[i]
}

save(pollen2k, file=paste0('data/pollen2k_', Sys.Date(), '.rdata'))
# 
# load('eda/rdata/pollen2k.rdata')
n = length(pollen2k)

# download the data
pollen_compiled = list()
for (i in 1:length(pollen2k)){  
   x = compile_list_neotoma(pollen2k[[i]][[1]], 'Stepps')
   pollen_compiled[[i]] = compile_list_stepps(x, list.name=list_name, pollen.equiv.stepps, cf = TRUE, type = TRUE)
   pollen_compiled[[i]]$counts_raw = pollen2k[[i]]$counts
   pollen_compiled[[i]]$dataset$state = pollen2k[[i]]$dataset$state
   pollen_compiled[[i]]$dataset$elev  = pollen2k[[i]]$dataset$elev
}

save(pollen_compiled, file=paste0('data/pollen_compiled', list_name, '_', Sys.Date(), '.rdata'))

#############################################################################
# Get the correct ages: choose NAPD and radiocarbon if present
#############################################################################

pollen_ts  = list()
age_models = vector(length=0)
types      = vector(length=0)

for (i in 1:length(pollen_compiled)){
  
  print(i)
  
  id       <- pollen_compiled[[i]]$dataset$dataset.meta$dataset.id
  lat      <- pollen_compiled[[i]]$dataset$site.data$lat
  long     <- pollen_compiled[[i]]$dataset$site.data$long
  elev     <- pollen_compiled[[i]]$dataset$site.data$elev
  state    <- pollen_compiled[[i]]$dataset$state
  sitename <- gsub("[ ']", "", as.vector(pollen_compiled[[i]]$dataset$site.data$site.name))
  
  # chronology name column indices
  idx         = which(colnames(pollen_compiled[[i]]$sample.meta) == "chronology.name")
  chrons_full = sapply(pollen_compiled[[i]]$sample.meta[1, idx], as.character)

  chrons_split = strsplit(chrons_full, split=" ")
  chron        = unname(sapply(chrons_split, "[[", 1))
  type         = as.character(pollen_compiled[[i]]$sample.meta[1, idx + 1])
  
  age_models = c(age_models, chron) 
  types = c(types, type)
  print(chron)

  ages = pollen_compiled[[i]]$sample.meta[, 'age']
  
  meta = list(id = id, sitename = sitename, lat = lat, long = long, state = state, elev = elev, model = chron)
  
  pollen_ts[[i]] = list(meta = meta, counts = data.frame(ages, pollen_compiled[[i]]$counts))
#   pollen_ts_nap[[i]] = list(meta = meta, counts = data.frame(ages, pollen_compiled[[i]]$counts_pls, pollen_compiled[[i]]$counts_nap))
}

#write it to a csv file: pollen_ts.csv
pollen_ts_df = list()
for (i in 1:length(pollen_ts)){
  print(i)
  counts = pollen_ts[[i]]$counts
  meta   = t(replicate(nrow(counts),unlist(pollen_ts[[i]]$meta)))
  pollen_ts_df = rbind(pollen_ts_df, cbind(meta, counts))
}
write.table(pollen_ts_df, file=paste('data/pollen_ts_', Sys.Date(), '.csv', sep=''), quote=FALSE, row.names=FALSE)



#############################################################################
# Make settlement era data set: ~120 ybp (=1830)
#############################################################################

build_pollen_se <- function(pollen_ts){
  
  n = length(pollen_ts)
  
  pollen_se     = list()
  nages         = vector(length=n)
  closest_age   = vector(length=n) # find the age closest to 120 ybp
  before120_age = vector(length=n) # find the age closest to but prior to 120 ypb
  
  for (i in 1:n){
    
    idx = which((pollen_ts[[i]]$counts$ages < (120+500)) & (pollen_ts[[i]]$counts$ages > (120-30)))
    nages[i] = length(idx)
    
#     pollen_se[[i]] = list(meta = pollen_ts[[i]]$meta, counts = pollen_ts[[i]]$counts[idx,])
    
    counts = pollen_ts[[i]]$counts[idx,]
    if (nrow(counts) > 0){
      meta   = t(replicate(nrow(counts),unlist(pollen_ts[[i]]$meta)))
      pollen_se = rbind(pollen_se, cbind(meta, counts))
    }
    
#     #do some checking
#     idx_closest_age = 0
#     closest_age[i] = NA
#     before120_age[i] = NA
#     if (length(pollen_se[[i]]$counts$ages) > 0){
#       idx_closest_age = which.min(abs(pollen_se[[i]]$counts$ages - 120))
#       closest_age[i] = pollen_se[[i]]$counts$ages[idx_closest_age]
#       
#       idx_before120    = which((pollen_se[[i]]$counts$ages - 120)>0)
#       
#       if (length(idx_before120)>0 ){
#         before120_age[i] = min(pollen_se[[i]]$counts$ages[idx_before120])
#       }
#     }
#     
#     print(pollen_se[[i]]$counts$ages)
#     print(closest_age[i])
  }
  
  return(pollen_se)
}

neo_se = build_pollen_se(pollen_ts)  
 
## CLH
clh    = read.csv('data/pollen_counts_calcote.csv', stringsAsFactors=FALSE)
clh_se = get_calcote(clh)

#####################################################################################
# build calibration data set from sites with only core-top and pre-settlement
#####################################################################################

clh.sites  <- read.csv('data/hotchkiss_lynch_calcote_meta_v0.1.csv', stringsAsFactors = FALSE, sep=',')
clh.counts <- read.csv('data/hotchkiss_lynch_calcote_counts_v0.1.csv', stringsAsFactors = FALSE, sep=',')

clh.sites$name <- gsub(" ","", clh.sites$name, fixed=TRUE)
clh.counts$name <- gsub(" ","", clh.counts$name, fixed=TRUE)

clh_pre = get_clhpre_se(clh.sites, clh.counts, list_name=list_name)$clh_pre

remove_factors <- function(df) {
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df
}

pollen_se_all = rbind(remove_factors(neo_se), remove_factors(clh_se), remove_factors(clh_pre))

write.table(pollen_se_all, file=paste('data/pollen_se_', Sys.Date(), '.csv', sep=''), quote=FALSE, row.names=FALSE)

# #############################################################################
# # dataset of counts summed over 120+500 to 120-30
# #############################################################################
# 
# pollen_se_sum = list()
# k = 0
# 
# for (i in 1:length(pollen_se)){
#   
#   meta   = pollen_se[[i]]$meta
#   counts_all = pollen_se[[i]]$counts[,-1]
#   
#   if (nrow(counts_all) != 0){
#     k = k + 1
#     counts = colSums(pollen_se[[i]]$counts[,-1])
#     pollen_se_sum[[k]] = list(meta=meta, counts=counts)
#   }
# }
# 
# 
# save(pollen_se_sum, file=paste('data/pollen/pollen_se_sum_', Sys.Date(), '.rdata', sep=''))
# 
# #write it to a csv file: pollen_ts.csv
# pollen_se_sum_df = list()
# for (i in 1:length(pollen_se_sum)){
#   print(i)
#   counts = pollen_se_sum[[i]]$counts
#   meta   = unlist(pollen_se_sum[[i]]$meta)
#   pollen_se_sum_df = rbind(pollen_se_sum_df, c(meta, counts))
#   
# }
# write.table(pollen_se_sum_df, file=paste('data/pollen/pollen_se_sum_', Sys.Date(), '.csv', sep=''), quote=FALSE, row.names=FALSE)