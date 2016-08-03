library(fields)

# pls       = read.csv('data/western_comp_v0.5.csv')
pls       = read.csv('data/composition/pls/western_comp_v0.6-2.csv')
pls.equiv = read.csv('pls.equiv.csv', stringsAsFactors=FALSE)

comp      = read.csv('data/composition/stat/composition_v0.5.csv')
centers_veg = comp[, c('x', 'y')]

list_name = 'must_have'

# path_data  = '../stepps-calibration/r/dump'
# suff_dat   = '12taxa_mid_comp_ALL_v0.3'
# # load composition data; r and centers_veg
# load(sprintf('%s/cal_data_%s.rdata', path_data, suff_dat))

keep_x = pls$x %in% centers_veg[,1]
keep_y = pls$y %in% centers_veg[,2]
keep_coord = keep_x & keep_y

pls = pls[keep_coord,]

dist_mat = rdist(as.matrix(centers_veg), as.matrix(pls[,1:2]))
dist_mat[dist_mat < 1e-6] = 0

idx = apply(dist_mat, 2, function(x) if (any(x == 0)) {which(x==0)}else{0})
pls2 = pls[idx != 0, ]

pls = pls2
# idx = idx[idx != 0]
# pls = pls[idx,]

# check
# plot(pls$x, pls$y)

# foo = pls[,3:ncol(pls)]

colnames(pls)[3:ncol(pls)] = pls.equiv[match(colnames(pls)[3:ncol(pls)], pls.equiv[,1]), list_name]

# bar = pls[,3:ncol(pls)]

sapply(unique(colnames(pls)[3:ncol(pls)]), function(x) print(x))


# res2 <-  sapply(unique(colnames(pls)[3:ncol(pls)]), function(x) 
#                rowMeans(pls[,colnames(pls)[3:ncol(pls)]== x] ))
# 
# res2 <-  sapply(unique(colnames(foo)), function(x) 
#   
#                  rowMeans(foo[,colnames(foo)== x]))
               

# aggregate(pls[,3:ncol(pls)], by=c(unique(colnames(pls)[3:ncol(pls)])), FUN=rowSums)


# sapply(unique(names(pls)[3:ncol(pls)]), 
#        function(x) rowSums( pls[ , grep(x, names(pls)[3:ncol(pls)]), drop=FALSE]) )

col_names = colnames(pls)[3:ncol(pls)]
nms = unique(col_names[!is.na(col_names)])

# sapply(nms, function(x) if (sum(col_names  == nms[1], na.rm=TRUE) > 1) {rowSums(pls[, x])} else {pls[,x]} )

sum_cols <- function(x){
  if (length(grep(x, names(pls))) > 1 ){
    rowSums(pls[,grep(x, names(pls))])
  } else {
    pls[,grep(x, names(pls))]
  }
}

pls_agg = sapply(nms, sum_cols)
pls_agg = pls_agg[, sort(colnames(pls_agg))]
pls = cbind(pls[,1:2], pls_agg)

# % of each taxon
colSums(pls[,3:ncol(pls)])/sum(colSums(pls[,3:ncol(pls)]))

write.table(pls, file='data/composition/pls/pls_umw_v0.6.csv', sep=',')
