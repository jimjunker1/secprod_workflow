## INS_JULIAN, converting and inserting JULIAN dates, Author: Jim Junker
## Last Updated: Sept. 2016
##
##Purpose: This code automates the plotting of length-frequency histograms of specified site, habitat, and taxon for a specified date range.
##	This function is meant to be used with a matrix of invertebrate size-classes(mm) and abundance (m^-2). There are five preceeding columns
##	SITE, SAMPLE, DATE, HABITAT, TAXON.

##	Date should maintain the following form with exact headings on columns 1-5

###	SITE	SAMPLE	DATE		HABITAT	TAXON		0.25	0.5	1	1.5
#	HVER	2		08/02/2011	COBBLE	OSTRACOD	0	0	2	2
#	HVER	2		08/02/2011	COBBLE	MIDGE		0	0	10	7
#	HVER	2		08/02/2011	COBBLE	OSTRACOD	0	0	6	10
#	HVER	2		08/02/2011	COBBLE	TUBIFICID	0	10	5	0
#	HVER	4		08/02/2011	COBBLE	LIMNOPHORA	0	0	0	0
#	.
#	.
#	.

## Notes:
##	SITE can handle an unlimited number of sites
##	SAMPLE number is detected and can handle unequal sample numbers across site/date/habitat/etc.
##	DATE This is the date of sampling. This must be in the form of MM/DD/YYY (sample date: "1/5/2008", "9/15/2007", etc.). This maintains the format for further
##		secondary production estimation.
##	HABITAT For habitat specific production. This is necessary for further secondary production estimation. Can have unlimited number of habitat types
##	TAXON This is the taxon information used for each group to summarize. This code will summarize the taxon with exactly matching names within samples.
##	Column 6 - n (as many as needed) These headings denote mm size classes. The values in these columns are abundances sample (m^2).
##		There is no requirement for these to be in order or in constant intervals. All columns should be filled in with counts or zeros. Truly missing data
##		should be filled with "NA". There is no need to include all non-occurrences of taxa, but will work fine if it does. 
##
##	When loading in the data file be sure to turn name checking off otherwise  
##	(i.e. DATA <- read.table(C:\DATA.txt, header = T, sep = "\t", quote = "", strip.white = T, row.names = F, check.names = F))
## 
ins_julian <- function(DATA, file.name,...) {
  
  ##Required packages 'chron'
  if(!require("pacman")) install.packages("pacman")
  library(pacman)
  package.list <- c("chron")
  p_load(char = package.list, install = T)
  rm("package.list")

#Extract years, months, and days, and re-code them as numeric variables:
year <- as.numeric(as.character(years(chron(dates = as.character(DATA$DATE)))))
month <- as.numeric(months(chron(dates = as.character(DATA$DATE))))
day <- as.numeric(days(chron(dates = as.character(DATA$DATE))))

#Combine year, month, and day into a single julian date variable starting with January 1, 2006 (i.e.- Jan 1, 2006 = day1):
JULIAN <- julian(month, day, year, origin=c(month = 12, day = 31, year = 2010))

#Insert the julian date variable into the dataframe:
DATA <- data.frame(DATA[,1:2], JULIAN, DATA[,3:(dim(DATA)[2])])

#sink(paste(file.name,"_julian.txt", sep = ""))
write.table(DATA, file = paste(file.name,"_julian.txt"), sep ="\t", row.names = F)
#sink()
}
