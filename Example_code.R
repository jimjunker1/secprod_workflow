##This is the example for testing the scripts in:
##frac_merge
##len_freq
##

#Set working directory
setwd("C:/Users/Jim/Documents/Projects/Iceland/Bug Samples/Secondary Production/Secondary Production R code suite")

#run the source of scripts for function
source("C:/Users/Jim/Documents/Projects/Iceland/Bug Samples/Secondary Production/Secondary Production R code suite/frac_merge/frac_merge_function.txt")
#source("C:/Users/Jim/Documents/Projects/Iceland/Bug Samples/Secondary Production/Secondary Production R code suite/len_freq/len_freq_plot.txt")

#load in the raw .csv file
#bugs <- read.csv(file = './Hver_samples.csv', T, check.names = F)

#load in raw .txt file

bugs <- read.table(file = './Hver_samples.txt', header = T, sep = "\t", quote = "", strip.white = T, check.names = F)

#merge the raw data by fraction

frac_merge(bugs, "Sample")

#Now work with the merged data

bugs_merge <- read.table(file = './Sample_merged.txt', header= T, sep = "\t", quote = "", strip.white = T, check.names = F, row.names = 1) #currently this has row names dammit!!

bugs_sum <- aggregate(bugs_merge[c(names(bugs_merge[6:(dim(bugs_merge)[2])]))], by = bugs_merge[c("SITE","DATE","HABITAT","TAXON")], FUN = median)

year <- as.numeric(as.character(years(chron(dates = as.character(bugs_sum$DATE)))))
month <- as.numeric(months(chron(dates = as.character(bugs_sum$DATE))))
day <- as.numeric(days(chron(dates = as.character(bugs_sum$DATE))))

JULIAN <- julian(month, day, year, origin = c(month = 12, day = 31, year = 2010))

bugs_sum = data.frame(bugs_sum[,1:2], JULIAN, bugs_sum[,3:dim(bugs_sum)[2]], check.names = F)


#dates = sort(unique(bugs_sum$JULIAN))



#Sample code to create the length frequency plot for a single taxon
hver_midge = bugs_sum[which(bugs_sum$TAXON == "Midge 1"),] #subset data to a single taxa


#Build histogram for a single date, single taxon
hv_midge_Aug = hver_midge[which(hver_midge$DATE == "08/02/2011"),] #subset data to single date

hv_midge_prop = prop.table(hv_midge_Aug[5:(dim(hv_midge_Aug)[2])]) #calculates relative frequency of mm size classes

hv_midge_data = gather(hv_midge_prop, size, rel.freq, 1:(dim(hv_midge_prop)[2]))

hv_hist <- ggplot(hv_midge_data, aes(x = as.numeric(size), y = rel.freq)) + 
	geom_bar(stat = "identity") +
	xlab("Size(mm)") +
	ylab("Relative Frequency");hv_hist
###

## For all dates 
hv_midge_prop2 = t(apply(hver_midge[6:(dim(hver_midge)[2])],1, prop.table)) #applies prop.table() function across all rows and transposes it back to same dimensions of col names
hv_midge_col = hver_midge[,1:5] #pulls the columns to add in later

hv_prop_table = data.frame(hv_midge_col, hv_midge_prop2, check.names = F) #makes data frame of site, date, hab., taxon and size class relative frequency

head(hv_prop_table)

#now writing a loop function to plot each date len_freq hist
loop = function(x,...){
for( i in x$DATE){

date_data = gather(x, size, rel.freq, 5:(dim(x)[2]))
hist <- ggplot(date_data, aes(x = as.numeric(size), y = rel.freq)) +
	geom_bar(stat = "identity") +
	xlab("Size (mm)") +
	ylab("Relative Frequency"); hist
dev.new() 
}
}


dates = as.character(sort(unique(hv_prop_table$JULIAN)))

pltList = list()
for(i in dates){
	data = hv_prop_table[which(hv_prop_table$JULIAN == i),];
	date_data = gather(data, size, rel.freq, 6:(dim(data)[2]));
	v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
	ymax = max(date_data$rel.freq) + 0.05
	xmax = max(date_data$size[which(date_data$rel.freq != 0)])

	pltList[[i]] <- ggplot(date_data, aes(x = as.numeric(size), y = rel.freq)) +
	geom_bar(stat = "identity", width = .99) +
	geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
	scale_y_continuous(limits = c(0,as.numeric(ymax)), expand = c(0,0)) +
	scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
	labs(x = "Size (mm)", y = "Relative Frequency") +
	ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == i)])))
#do.call(plot_grid,pltList)
args.list = list(grobs = pltList, top = paste(date_data$TAXON))
#do.call(grid.arrange, pltList)
do.call(grid.arrange, args.list)
#print(multiplot(plotlist = pltList, ncol = 3))
print(i)
}  #this works, need to remove x lab and y labs and fit it better into page.

