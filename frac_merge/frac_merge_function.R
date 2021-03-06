## FRAC_MERGE, fraction consolidation code, Author: Jim Junker
## Last Updated: Sept. 2016
## code runtime for a data set with a single stream is currently 49 secs
##Purpose: This code automates consolidation of fractioned invertebrate samples for secondary production estimation. 
##	It sums invert size classes across small and large fractions within sites, dates, repilcates, habitats, and taxa while accounting for variable subsampling 
##	practices among split fractions. Further, it scales samples to a common areal basis from raw count-size class input data while being flexible for potential 
##	changing sampler area. It outputs a new .csv file with all sample's fractions summed and merged and on a areal basis.
##
##	Data should maintain the following form with exact column headings in columns 1-8
##
###	SITE	SAMPLE	FRACTION	SUBSAMPLE	AREA_COR	DATE		HABITAT	TAXON		0.25	0.5	1	1.5...
##	HVER	2		>1MM		1		0.09		08/02/2011	COBBLE	OSTRACOD	0	0	2	2
##	HVER	2		>1MM		1		0.09		08/02/2011	COBBLE	MIDGE		0	0	10	7
##	HVER	2		<1MM		16		0.09		08/02/2011	COBBLE	OSTRACOD	0	0	6	10
##	HVER	2		<1MM		16		0.09		08/02/2011	COBBLE	TUBIFICID	0	10	5	0
##	HVER	4		>1MM		1		0.09		08/02/2011	COBBLE	LIMNOPHORA	0	0	0	0
##	.
##	.
##	.
## Notes: 
##	SITE can handle an unlimited number of sites
##	SAMPLE number is detected and can handle unequal sample numbers across site/date/habitat/etc.
##	FRACTION only handles 2 fractions labled ">1MM" and "<1MM". A warning is given if a sample does not contain paired sample fractions, so empty samples must be input
##	SUBSAMPLE This represents the fraction. Counts will be multiplied by this number
##	AREA_COR This represents the proportion of 1 meter^2 covered with the sampler. Counts are divided by this proportion to scale to ind. * m^-2 
##	DATE This is the date of sampling. This must be in the form of MM/DD/YYY (sample date: "1/5/2008", "9/15/2007", etc.). This maintains the format for further
##		secondary production estimation.
##	HABITAT For habitat specific production. This is necessary for further secondary production estimation.
##	TAXON This is the taxon information used for each group to merge. This code will merge all taxon with exactly matching names within samples.
##	Column 9 - n (as many as needed) These headings denote mm size classes. The values in these columns are abundances in the sample/subsample.
##		There is no requirement for these to be in order or in constant intervals. All columns should be filled in with counts or zeros. Truly missing data
##		should be filled with "NA". There is no need to include all non-occurrences of taxa, but will work fine if it does.. 
##
##	When loading in the data file be sure to turn name checking off otherwise  
##	(i.e. DATA <- read.table(C:\DATA.txt, header = T, sep = "\t", quote = "", strip.white = T, check.names = F))
## 
##Required packages 'plyr', 'dplyr', 'tidyr', 'reshape2', 'chron'
if(!require("pacman")) install.packages("pacman")
library(pacman)
package.list <- c( "plyr", "dplyr", "tidyr", "reshape2", "chron")
p_load(char = package.list, install = T)
rm(package.list)
##
frac_merge <- function(DATA, file.name,...) {
	stm <- proc.time()
##
	DATA <- DATA[!is.na(apply(DATA[,9:(dim(DATA)[2])],1,sum)),]			#Removes NA's (missing data)
	row.names(DATA) <- NULL									#Resets row names to be sequential
##
	DATA_LONG <- gather(DATA, size_class, count, 9:(dim(DATA)[2]))		#Transforms the data set from wide to long
##
##Make sure there is a large and small fraction for each size, sample, date, habitat
	DATA_CHECK <- ddply(DATA_LONG, c("SITE", "SAMPLE", "DATE", "FRACTION"), summarize, length(FRACTION))	#get
#	DATA_CHECK <- ddply(DATA_CHECK, c("SITE", "SAMPLE", "DATE", "FRACTION"), summarize, LENGTH = length(FRACTION))
##
	FRAC_CHECK <- dcast(DATA_CHECK, SITE + SAMPLE + DATE ~ FRACTION, value.var = "FRACTION")
	colnames(FRAC_CHECK) <- c("SITE", "SAMPLE", "DATE", "SMALL", "LARGE")
	for(i in 1:nrow(FRAC_CHECK)){
		if(is.na(FRAC_CHECK$SMALL[[i]] != "<1mm" | FRAC_CHECK$LARGE[[i]] != ">1mm")){
	  		print("ERROR: there are samples unpaired sample fractions")
	  		print(rownames(FRAC_CHECK[i,]))
	}}
##Continue with the fraction merging
	DATA_LONG.sub <- ddply(DATA_LONG, c("SITE", "SAMPLE", "FRACTION", "AREA_COR", "DATE", "HABITAT", "TAXON", "size_class"), transform, count = count*SUBSAMPLE)#Applies a subsample correction for each sample and fraction
	DATA_LONG.area <- ddply(DATA_LONG.sub, c("SITE", "SAMPLE", "DATE", "HABITAT", "TAXON", "size_class"), transform, count = count/AREA_COR)				#Applies areal correction for each sample
	DATA_LONG.sum <- ddply(DATA_LONG.area, c("SITE", "SAMPLE", "DATE", "HABITAT", "TAXON", "size_class"), summarize, count_fin = sum(count))	#Sums fractions for each size class to get to density
	DATA_FIN <- spread(DATA_LONG.sum, size_class, count_fin)
	#sink(paste(file.name,"_merged.txt", sep = ""))
	write.table(DATA_FIN, file = paste(file.name,"_merged.txt"), sep = "")
	#sink()
	total.time <- proc.time() - stm
	print(paste("Total run time = ", total.time[[1]]))
}
