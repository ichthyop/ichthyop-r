compute_recruitment <- function(){

	library(ncdf)

	# The directory that contains the series of netcdf input files
	dirpath <- 'F:\\SAfE_SOUTHCOAST\\IchthyopOutputs\\Ichthyop2.1\\Child\\Serial_transport_1\\'
	
	# In case one wishes to consider only a subset of all drifters
	# Index of the first and last drifters
	firstdrifter <- 1
	lastdrifter <- 10000

	# The time record at which to compute recruitment
	computeattime <- 13

	# The number of release zones
	nbreleasezones <- 5

	# The index of the recruitment zone for which recruitment is computed
	recruitmentzone <- 1

	# The maximum value of the transport success plot (in %)
	ymax <- 50

	# An inner function that computes recruitment for one file
	compute_recruitment_file <- function(filename){

		# An inner function that computes year and day from time in seconds
		compute_yearday <- function(time){
			nbdays = 1+time/86400
			year = 1+as.integer(nbdays/360)
			day = as.integer(nbdays-360*(year-1))
			if (day < 10) day <- paste('0',day,sep='')
			if (as.numeric(day) < 100) day <- paste('0',day,sep='')
			return(c(year,day))
		}

		nc <- open.ncdf(filename)

		# Gets the value of time of release
		t0 <- get.var.ncdf(nc,'time',1,1)
		# Computes year and day of release
		yearday <- compute_yearday(t0)
		# Gets the values of min and max depth of release and concatenates them
		mindepth <- att.get.ncdf(nc,0,'release_depth_min')$value
		maxdepth <- att.get.ncdf(nc,0,'release_depth_max')$value
		depth <- paste(mindepth,maxdepth,sep='-')
		# Gets the value of replica number
		replica <- att.get.ncdf(nc,0,'replica')$value

		nbdrifter <- lastdrifter-firstdrifter+1
		# Gets the value of recruited for the recruitment zone considered for all drifters at time of computation
		recruited <- get.var.ncdf(nc,'recruited',c(recruitmentzone,firstdrifter,computeattime),c(1,nbdrifter,1))
		# Gets the value of release zone for all drifters
		releasezone <- get.var.ncdf(nc,'zone',c(firstdrifter,1),c(nbdrifter,1))
		# Calculates the number of recruits from every release zone
		recruitnb <- hist(recruited*releasezone,seq(0,nbreleasezones+1)-0.5,plot=FALSE)$counts[2:(nbreleasezones+1)]
		# Calculates the number of released from every release zone
		releasenb <- hist(releasezone,seq(0,nbreleasezones+1)-0.5,plot=FALSE)$counts[2:(nbreleasezones+1)]
		close.ncdf(nc)
		# returns a collage of columns, i.e., a table, that looks like the following
		# releasenb1 recruitnb1 1 year day depth replica
		# releasenb2 recruitnb2 2 year day depth replica
		# releasenb3 recruitnb3 3 year day depth replica
		# ...
		return(cbind(releasenb,recruitnb,seq(1,nbreleasezones),rep(yearday[1],nbreleasezones),
			rep(yearday[2],nbreleasezones),rep(depth,nbreleasezones),rep(replica,nbreleasezones)))
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

	# Computes the %age of recruitment for every release area, year, day, depth and replica
	dataarea <- 100*tapply(as.numeric(dataset[,2]),dataset[,3],sum)/tapply(as.numeric(dataset[,1]),dataset[,3],sum)
	datayear <- 100*tapply(as.numeric(dataset[,2]),dataset[,4],sum)/tapply(as.numeric(dataset[,1]),dataset[,4],sum)
	dataday <- 100*tapply(as.numeric(dataset[,2]),dataset[,5],sum)/tapply(as.numeric(dataset[,1]),dataset[,5],sum)
	datadepth <- 100*tapply(as.numeric(dataset[,2]),dataset[,6],sum)/tapply(as.numeric(dataset[,1]),dataset[,6],sum)
	datareplica <- 100*tapply(as.numeric(dataset[,2]),dataset[,7],sum)/tapply(as.numeric(dataset[,1]),dataset[,7],sum)

	par(mfrow=c(4,2))

	# Makes the corresponding plots
	areaplot <- barplot(dataarea,xlab='Release area',ylab='Transport success (%)',ylim = c(0,ymax))
	#areaplot <- barplot(dataarea,xlab='Release area',ylab='Transport success (%)',ylim = c(0,ymax),names.arg = c("EABin","EABoff","CABin","CABoff","WAB"))
	dayplot <- barplot(dataday,xlab='Release day',ylab='Transport success (%)',ylim = c(0,ymax))
	#dayplot <- barplot(dataday,xlab='Release month',ylab='Transport success (%)',ylim = c(0,ymax),names.arg = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
	yearplot <- barplot(datayear,xlab='Release year',ylab='Transport success (%)',ylim = c(0,ymax))
	depthplot <- barplot(datadepth,xlab='Release depth (m)',ylab='Transport success (%)',ylim = c(0,ymax))
	replicateplot <- barplot(datareplica ,xlab='Replicate',ylab='Transport success (%)',ylim = c(0,ymax))

	#recruitprop <- 100*as.numeric(dataset[,2])/as.numeric(dataset[,1])
	#dataset <- as.data.frame(dataset)
	#colnames(dataset) <- c('NumberReleased','NumberRecruited','ReleaseArea','Year','Day','Depth','Replicate')
	#mod <- lm(recruitprop ~ factor(ReleaseArea) + factor(Day) + factor(Year) + factor(Depth) + factor(Replicate)
	#			+ factor(ReleaseArea):factor(Day) + factor(ReleaseArea):factor(Year) + factor(ReleaseArea):factor(Depth)
	#			+ factor(Day):factor(Year) + factor(Day):factor(Depth) + factor(Year):factor(Depth), data = dataset)
	#aov <- anova(mod)
	#print(aov)
	#print(100 * aov[2] / sum(aov[2]))
	#interaction.plot(dataset$Day,dataset$ReleaseArea,recruitprop,fixed=TRUE,xlab='Release day',ylab='Transport success (%)',lty=1,col=seq(1,length(areaplot)))
	#interaction.plot(dataset$Depth,dataset$ReleaseArea,recruitprop,fixed=TRUE,xlab='Release depth (m)',ylab='Transport success (%)',lty=1,col=seq(1,length(areaplot)))
	#interaction.plot(dataset$Day,dataset$Year,recruitprop,fixed=TRUE,xlab='Release day',ylab='Transport success (%)',lty=1,col=seq(1,length(yearplot)))
}
