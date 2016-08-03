#install “vegan” and “rioja” LATEST packages (rioja must be 0.9-6 or later)
#read in data
Devil <- read.csv("DevilPollen.csv")
DevilPrct <- Devil #calculate percentages
DevilPrct[, 7:32] <- Devil[, 7:32]/(Devil[, 6])*100 #change column inputs
#'7' above is the row species start at in the .csv file
#in DevilPollen there are 32 species variables
#countallgrains in csv are total grains (not all necessarily used in your pollen sum)
#countgrainsshownhere in csv in column 6 is your pollen sum used to calculate percentages

# construct pollen diagram
spec <- DevilPrct [, 7:32] # spec for species pollen percentages
strat.plot (spec, yvar = DevilPrct$Age, exag=TRUE, scale.percent = TRUE, y.rev = TRUE, ylabel = "Age (cal years BP)", plot.bar = FALSE, plot.poly = TRUE)

#can use age or depth for y axis, just need to call that row with $

# calculate distance 
diss <- dist (sqrt(spec/100))

# add clustering
clust <- chclust (diss, method = "coniss")
bstick (clust)

#add zones
x <- strat.plot (spec, scale.percent = TRUE, graph.widths=0, yvar = DevilPrct$Age, y.rev = TRUE, ylabel = "Age (cal years BP)", plot.bar = FALSE, plot.poly = TRUE, 
                 clust = clust, cex.xlabel=0.8, srt.xlabel=65, exag=TRUE)
addClustZone (x, clust, 4, col = "red")

# adjust the diagram position and size, and export the diagram

x <- strat.plot (spec, xLeft = 0.09, xRight = 0.98, yvar = DevilPrct$Age, y.rev = TRUE, ylabel = "Age (cal years BP)", plot.bar = FALSE, plot.poly = TRUE, 
                 scale.percent = TRUE, clust = clust, cex.xlabel=0.7, srt.xlabel=65, exag=TRUE)
addClustZone (x, clust, 4, col = "red")


#for more ref only!!
help(strat.plot)
