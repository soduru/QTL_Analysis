# R script for analyzing QTL Sugiyama data
# Srijan Oduru
# sugiyamashort.qtl script
# October 10 2021
#
# clean things up
rm(list=ls())
setwd("~/Desktop/R Folder")
#
# load the QTL library
# NOTE!  I first had to INSTALL the library using: install.packages("qtl")
# Now I can use the package qtl
library(qtl)
#
# Load the data!
sugiyama <- read.cross("csv", file="sugiyamashort.csv", genotypes=c("C", "H", "B"), 
                    na.strings="-", alleles=c("B", "C"))
# Sometimes the genetic markers are too close. Jittermap will move them apart slightly so my results are better.
jittermap(sugiyama)
#  A summary of the cross gives me some basic data.  Nice!
summary(sugiyama)
# I need to see what phenotypes are in the dataset.  names does that for me
names(sugiyama)
# take a look at my data, make sure it's pretty clean
#  I should NOT see any really big red spots, especially in the bottom right corner under the diagonal
sugiyama <- est.rf(sugiyama)
plotRF(sugiyama)
# It's nice to see my genetic map -- all of the horizontal lines are genetic markers that have been inserted.
plot.map(sugiyama)
#  It's often the case that I have missing data -- plot.missing shows me where it is
plotMissing(sugiyama)
# a histogram of my BP phenotype -- if I get a "normal" distribution (a bell-shaped curve), then I can be pretty confident of the data
hist(sugiyama$pheno$BP_final, main="Histogram of BP")
hist(sugiyama$pheno$HR_final, main="Histogram of Resting HR")
hist(sugiyama$pheno$heart_wt, main="Histogram of Heart Weight")
#  I'm lazy....I hate to type "cross$pheno$bp every time I need to study the BP, so I give it a short name.



#
# Now I'm going to generate a mainscan.  First, I calculate what the scan SHOULD look like, so I'm going to calculate a genetic probability map.
sugiyama <- calc.genoprob(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane",
                       stepwidth="fixed")
#
# Run a simulated geno probability calculations
sugiyama <- sim.geno(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane", 
                  stepwidth="fixed")
# Perform the mainscan for the QTL
sugiyama.scanBP <- scanone(sugiyama, pheno.col=3, model="normal", method="em")
#  I'm only going to run this for 100 "permulations" -- typically you do 500-1000, but that takes a LONG time.
sugiyama.scanBP.perm <- scanone(sugiyama, pheno.col=3, model="normal", method="em", n.perm=100)

#
# plot the mainscan
plot(sugiyama.scanBP, main="Mainscan of BP")
#  I'm putting threshold likes at 63% confidence, 90% confidence, and 95% confidence.
thresh <- summary(sugiyama.scanBP.perm, alpha=c(0.63, 0.10, 0.05))
abline(h=thresh[1], col="red")
abline(h=thresh[2], col="blue")
abline(h=thresh[3], col="green")

#2nd iteration of mainscan generation and plotting resting heart rate
sugiyama <- calc.genoprob(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane",
                          stepwidth="fixed")
sugiyama <- sim.geno(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane", 
                     stepwidth="fixed")

sugiyama.scanHR_final <- scanone(sugiyama, pheno.col=4, model="normal", method="em")
sugiyama.scanHR_final.perm <- scanone(sugiyama, pheno.col=4, model="normal", method="em", n.perm=100)

plot(sugiyama.scanHR_final, main="Mainscan of HR")

thresh <- summary(sugiyama.scanHR_final.perm, alpha=c(0.63, 0.10, 0.05))
abline(h=thresh[1], col="red")
abline(h=thresh[2], col="blue")
abline(h=thresh[3], col="green")
#
#3rd iteration of mainscan generation and plotting heart weight
sugiyama <- calc.genoprob(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane",
                          stepwidth="fixed")
sugiyama <- sim.geno(sugiyama, step=2.0, off.end=0.0, error.prob=1.0e-4, map.function="haldane", 
                     stepwidth="fixed")

sugiyama.scanheart_wt <- scanone(sugiyama, pheno.col=6, model="normal", method="em")
sugiyama.scanheart_wt.perm <- scanone(sugiyama, pheno.col=6, model="normal", method="em", n.perm=100)

plot(sugiyama.scanheart_wt, main="Mainscan of Heart Weight")

thresh <- summary(sugiyama.scanheart_wt.perm, alpha=c(0.63, 0.10, 0.05))
abline(h=thresh[1], col="red")
abline(h=thresh[2], col="blue")
abline(h=thresh[3], col="green")
#
#  I'd like to see a text-based output of my scan
summary(sugiyama.scanBP, perm=sugiyama.scanBP.perm, lodcolumn=1, alpha=0.05)
summary(sugiyama.scanHR_final, perm=sugiyama.scanHR_final.perm, lodcolumn=1, alpha=0.05)
summary(sugiyama.scanheart_wt, perm=sugiyama.scanheart_wt.perm, lodcolumn=1, alpha=0.05)
#
# do an effect plot
# once you see an effect plot, you'll understand what it does!
# first effect plot (bp)
mname1 <- find.marker(sugiyama, chr=7, pos=48.7)
effectplot(sugiyama, pheno.col=3, mname1=mname1)
# second effect plot (bp)
mname2 <- find.marker(sugiyama, chr=15, pos=12.0)
effectplot(sugiyama, pheno.col=3, mname1=mname2)
# third effect plot (hr)
mname3 <- find.marker(sugiyama, chr=2, pos=59.8)
effectplot(sugiyama, pheno.col=4, mname1=mname3)





# All done
#EOF