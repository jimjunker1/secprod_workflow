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


#This code needs to be integrated into the FRAC_MERGE script.
#Or more importantly should be made into its own function to adjust dates to Julian
#across all data.frames used to production analysis
year <- as.numeric(as.character(years(chron(dates = as.character(bugs_sum$DATE)))))
month <- as.numeric(months(chron(dates = as.character(bugs_sum$DATE))))
day <- as.numeric(days(chron(dates = as.character(bugs_sum$DATE))))

JULIAN <- julian(month, day, year, origin = c(month = 12, day = 31, year = 2010))

bugs_sum = data.frame(bugs_sum[,1:2], JULIAN, bugs_sum[,3:dim(bugs_sum)[2]], check.names = F)


#dates = sort(unique(bugs_sum$JULIAN))
######


#Sample code to create the length frequency plot for a single taxon
hver_midge = bugs_sum[which(bugs_sum$TAXON == "Midge 1"),] #subset data to a single taxa


#Build histogram for a single date, single taxon. don't run this because it doesn't work after 
#Julian date entry
#hv_midge_Aug = hver_midge[which(hver_midge$DATE == "08/02/2011"),] #subset data to single date

#hv_midge_prop = prop.table(hv_midge_Aug[5:(dim(hv_midge_Aug)[2])]) #calculates relative frequency of mm size classes

#hv_midge_data = gather(hv_midge_prop, size, rel.freq, 1:(dim(hv_midge_prop)[2]))

#hv_hist <- ggplot(hv_midge_data, aes(x = as.numeric(size), y = rel.freq)) + 
#	geom_bar(stat = "identity") +
#	xlab("Size(mm)") +
#	ylab("Relative Frequency");hv_hist
###

## For all dates 
hv_midge_prop2 = t(apply(hver_midge[6:(dim(hver_midge)[2])],1, prop.table)) #applies prop.table() function across all rows and transposes it back to same dimensions of col names
hv_midge_col = hver_midge[,1:5] #pulls the columns to add in later

hv_prop_table = data.frame(hv_midge_col, hv_midge_prop2, check.names = F) #makes data frame of site, date, hab., taxon and size class relative frequency

head(hv_prop_table)


#####Another try to plot all dates on single plot. Need to save file with 
#species names rather than making it 
dates = as.character(sort(unique(hv_prop_table$JULIAN)))

pltList = list()
for(i in dates) {
  stm <- proc.time()
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
	theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
	ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == i)])))
#do.call(plot_grid,pltList)
args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = "Relative Frequency", bottom = " Size (mm)")
#do.call(grid.arrange, pltList)
do.call(grid.arrange, args.list)
#print(multiplot(plotlist = pltList, ncol = 3))
#print(i)
total.time <- proc.time() - stm
print(paste("Run time=", round(total.time[[1]],digits = 2)))
}  #this works, need to remove x lab and y labs and fit it better into page.

#######

# Ideally,  I need to pu this into a function that automatically makes these plots so 
# I can use it multiple times in len_freq function
#here is the addition of a second taxon
hver_midge2 = bugs_sum[which(bugs_sum$TAXON == "Midge 1" | bugs_sum$TAXON == "Snail"),] #subset data to a single taxa

## For all dates 
hv_midge_prop3 = t(apply(hver_midge2[6:(dim(hver_midge2)[2])],1, prop.table)) #applies prop.table() function across all rows and transposes it back to same dimensions of col names
hv_midge_col2 = hver_midge2[,1:5] #pulls the columns to add in later

hv_prop_table2 = data.frame(hv_midge_col2, hv_midge_prop3, check.names = F) #makes data frame of site, date, hab., taxon and size class relative frequency
hv_prop_table2[is.na(hv_prop_table2)] <- 0
head(hv_prop_table2)

taxon = unique(levels(hv_prop_table2$TAXON))
for(k in taxon) {
  data1 = hv_prop_table2[which(hv_prop_table2$TAXON == k),]
  
for(i in dates) {
  stm <- proc.time()
  data = data1[which(data1$JULIAN == i),];
  date_data = gather(data, size, rel.freq, 6:(dim(data)[2]));
  v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
  ymax = max(date_data$rel.freq) + 0.05
  xmax = max(date_data$size[which(date_data$rel.freq != 0)])
  
  pltList[[i]] <- ggplot(date_data, aes(x = as.numeric(size), y = rel.freq)) +
    geom_bar(stat = "identity", width = .99) +
    geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
    scale_y_continuous(limits = c(0,as.numeric(ymax)), expand = c(0,0)) +
    scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    #labs(x = "Size (mm)", y = "Relative Frequency") +
    ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == i)])))
 
}  
  
   #do.call(plot_grid,pltList)
args.list = list(grobs = pltList, top = paste(date_data$TAXON), left = paste("Relative Frequency"), bottom = paste(" Size (mm)"))
  #do.call(grid.arrange, pltList)
  #plots <- do.call(grid.arrange, args.list)
  #print(multiplot(plotlist = pltList, ncol = 3))
  #print(i)
  total.time <- proc.time() - stm
  print(paste("Run time=", round(total.time[[1]],digits = 2)))
  pdf(file = paste(data$SITE[i],"_",data$Taxon[k],"_","plots.pdf"))
  #pdf("Plot.pdf")
  do.call(grid.arrange, pltList)
  #plots
  dev.off()
}

##########This below is for future. Forget it for now and just get shit done.######
##running in parallel to see if it speeds it up. Need package "foreach", "doParallel", "iterators"
cores = detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)


foreach(i= 1:dates) %do% {  #to run this parallel swith to %dopar%
  stm <- proc.time()
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
  total.time <- proc.time() - stm
  print(paste("Run time=", total.time[[1]]))
} 


#tryin something out here. Making the plotting a function so can run parallel or 
#just put in the len_freq function in multiple places
plot_loop = function(DATA,...) {
  stm <- proc.time()
  data = DATA[which(DATA$JULIAN == i),];
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
  total.time <- proc.time() - stm
  print(paste("Run time=", round(total.time[[1]],digits = 2)))
} 


foreach(i in dates) %do% plot_loop(i)

#old code to keep just in case. delete when done.#

#now writing a loop function to plot each date len_freq hist. this does not work for all dates.
loop = function(x,...){
  for( i in x$JULIAN){
    
    date_data = gather(x, size, rel.freq, 5:(dim(x)[2]))
    v = date_data$size[which(date_data$rel.freq == max(date_data$rel.freq))]
    ymax = max(date_data$rel.freq) + 0.05
    xmax = max(date_data$size[which(date_data$rel.freq != 0)])
    hist <- ggplot(date_data, aes(x = as.numeric(size), y = rel.freq)) +
      geom_bar(stat = "identity", width = 0.99) +
      geom_vline(xintercept = as.numeric(v), colour = "red", size = 1, linetype = "dashed") +
      scale_y_continuous(limits = c(0,as.numeric(ymax)), expand = c(0,0)) +
      scale_x_continuous(limits = c(0, as.numeric(xmax) + 1), expand = c(0,0)) +
      labs(x = "Size (mm)", y = "Relative Frequency") +
      ggtitle(paste(as.character(date_data$DATE[which(date_data$JULIAN == i)]))); hist
    dev.new() 
  }
}
###
