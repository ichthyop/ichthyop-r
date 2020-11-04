compute_recruitment_Ichthyop_v3 <- function(){

	library(ncdf)
	library(XML)

	# The directory that contains the series of netcdf input files
	dirpath <- 'C:\\Users\\Christophe\\IRD\\Ichthyop\\Ichthyop_stable-3\\stable-3\\stable-3\\output\\'
	
	# In case one wishes to consider only a subset of all drifters
	# Index of the first and last drifters
	firstdrifter <- 1
	lastdrifter <- 10000

	# The time record at which to compute recruitment
	computeattime <- 30

	# The number of release zones
	nbreleasezones <- 8

	# The index of the recruitment zone for which recruitment is computed
	recruitmentzone <- 1

	# The maximum value of the transport success plot (in %)
	ymax <- 6.5

	# An inner function that computes recruitment for one file
	compute_recruitment_file <- function(filename){

		# An inner function that computes year and day from time in seconds
		compute_yearday <- function(time){
			nbdays = 1+time/86400
			year = 1+as.integer(nbdays/360)
			day = as.integer(nbdays-360*(year-1))
			#if (day < 10) day <- paste('0',day,sep='')
			#if (as.numeric(day) < 100) day <- paste('0',day,sep='')
			return(c(year,day))
		}

		nc <- open.ncdf(filename)

		# Gets the value of time of release
		t0 <- get.var.ncdf(nc,'time',1,1)
		# Computes year and day of release
		yearday <- compute_yearday(t0)
		#Reads zonefile, gets the values of min and max depth of release and concatenates them
		#Warning all release zones should have the same min and max depths as the script returns the values found for the first zone only
		filezone <- att.get.ncdf(nc,0,'release.zone.zone_file')$value
		mindepth <- xmlValue(xmlRoot(xmlTreeParse(filezone))[["zone"]][["thickness"]][["upper_depth"]])
		maxdepth <- xmlValue(xmlRoot(xmlTreeParse(filezone))[["zone"]][["thickness"]][["lower_depth"]])
		depth <- paste(mindepth,maxdepth,sep='-')

		nbdrifter <- lastdrifter-firstdrifter
		# Gets the value of recruited for the recruitment zone considered for all drifters at time of computation
		recruited <- get.var.ncdf(nc,'recruited',c(recruitmentzone,firstdrifter,computeattime),c(1,nbdrifter,1))
		# Gets the value of release zone for all drifters
		releasezone <- get.var.ncdf(nc,'zone',c(1,firstdrifter,1),c(1,nbdrifter,1))
		releasezone <- releasezone +1
		# Calculates the number of recruits from every release zone
		recruitnb <- hist(recruited*releasezone,seq(0,nbreleasezones+1)-0.5,plot=FALSE)$counts[2:(nbreleasezones+1)]
		# Calculates the number of released from every release zone
		releasenb <- hist(releasezone,seq(0,nbreleasezones+1)-0.5,plot=FALSE)$counts[2:(nbreleasezones+1)]
		close.ncdf(nc)
		# returns a collage of columns, i.e., a table, that looks like the following
		# releasenb1 recruitnb1 1 year day depth
		# releasenb2 recruitnb2 2 year day depth
		# releasenb3 recruitnb3 3 year day depth
		# ...
		return(cbind(releasenb,recruitnb,seq(1,nbreleasezones),rep(yearday[1],nbreleasezones),
			rep(yearday[2],nbreleasezones),rep(depth,nbreleasezones)))
	}


	# An inner function that computes statistics of recruitment for a given factor
	compute_recruitment_stats <- function(released,recruited,factor){
		mean <- tapply(recruited,factor,sum)/tapply(released,factor,sum)
		var <- tapply(recruited^2,factor,sum)/tapply(released^2,factor,sum)-mean^2
		sem <- sqrt(var/table(factor))
		return(cbind(mean,sem))
	}


	# Gets filenames of all files in the dirpath directory 
	filenames <- list.files(path = dirpath, full.names = TRUE)
	# Computes recruitment for the first file 
	dataset <- compute_recruitment_file(filenames[1])
	# Computes recruitment for all subsequent files 
	if (length(filenames) >1){
		for (i in seq(2, length(filenames))){
			# Shows name of opened file on the console
			print(filenames[i])
			flush.console()
			# Computes recruitment data for file i
			data <- compute_recruitment_file(filenames[i])
			# Adds recruitment data computed for file i to those computed from all previous files
			dataset <- rbind(dataset,data)
		}
	}

	#print(dataset)

	# Computes stats (mean, std error of the mean) of recruitment for every release area, year, day, and depth
	
	dataarea_stats=compute_recruitment_stats(as.numeric(dataset[,1]),as.numeric(dataset[,2]),as.numeric(dataset[,3]))
	dataarea_mean=dataarea_stats[,1]
	dataarea_sem=dataarea_stats[,2]
	
	datayear_stats=compute_recruitment_stats(as.numeric(dataset[,1]),as.numeric(dataset[,2]),as.numeric(dataset[,4]))
	datayear_mean=datayear_stats[,1]
	datayear_sem=datayear_stats[,2]

	dataday_stats=compute_recruitment_stats(as.numeric(dataset[,1]),as.numeric(dataset[,2]),as.numeric(dataset[,5]))
	dataday_mean=dataday_stats[,1]
	dataday_sem=dataday_stats[,2]

	datadepth_stats=compute_recruitment_stats(as.numeric(dataset[,1]),as.numeric(dataset[,2]),dataset[,6])
	datadepth_mean=datadepth_stats[,1]
	datadepth_sem=datadepth_stats[,2]

	# Makes the corresponding plots
	par(mfrow=c(2,2))

	areaplot <- barplot(100*dataarea_mean,xlab='Release area',ylab='Transport success (%)',ylim = c(0,ymax))
	arrows(areaplot,100*(dataarea_mean+dataarea_sem),areaplot,100*(dataarea_mean-dataarea_sem),angle=90,code=3,length=0.05)
	dayplot <- barplot(100*dataday_mean,xlab='Release day',ylab='Transport success (%)',ylim = c(0,ymax))
	arrows(dayplot,100*(dataday_mean+dataday_sem),dayplot,100*(dataday_mean-dataday_sem),angle=90,code=3,length=0.05)
	yearplot <- barplot(100*datayear_mean,xlab='Release year',ylab='Transport success (%)',ylim = c(0,ymax))
	arrows(yearplot,100*(datayear_mean+datayear_sem),yearplot,100*(datayear_mean-datayear_sem),angle=90,code=3,length=0.05)
	depthplot <- barplot(100*datadepth_mean,xlab='Release depth (m)',ylab='Transport success (%)',ylim = c(0,ymax))
	arrows(depthplot,100*(datadepth_mean+datadepth_sem),depthplot,100*(datadepth_mean-datadepth_sem),angle=90,code=3,length=0.05)

	#recruitprop <- 100*as.numeric(dataset[,2])/as.numeric(dataset[,1])
	#dataset <- as.data.frame(dataset)
	#colnames(dataset) <- c('NumberReleased','NumberRecruited','ReleaseArea','Year','Day','Depth')
	#mod <- lm(recruitprop ~ factor(ReleaseArea) + factor(Day) + factor(Year) + factor(Depth)
	#			+ factor(ReleaseArea):factor(Day) + factor(ReleaseArea):factor(Year) + factor(ReleaseArea):factor(Depth)
	#			+ factor(Day):factor(Year) + factor(Day):factor(Depth) + factor(Year):factor(Depth), data = dataset)
	#aov <- anova(mod)
	#print(aov)
	#print(100 * aov[2] / sum(aov[2]))
	#interaction.plot(dataset$Day,dataset$ReleaseArea,recruitprop,fixed=TRUE,xlab='Release day',ylab='Transport success (%)',lty=1,col=seq(1,length(areaplot)))
	#interaction.plot(dataset$Depth,dataset$ReleaseArea,recruitprop,fixed=TRUE,xlab='Release depth (m)',ylab='Transport success (%)',lty=1,col=seq(1,length(areaplot)))
	#interaction.plot(dataset$Day,dataset$Year,recruitprop,fixed=TRUE,xlab='Release day',ylab='Transport success (%)',lty=1,col=seq(1,length(yearplot)))
}