plot_dens <- function(){

	library(ncdf)

	# The netcdf input file
	ncfile <- 'F:\\SAfE_SOUTHCOAST\\IchthyopOutputs\\serial_simu_0000.nc'

	# In case one wishes to consider only a subset of all drifters
	# Index of the first and last drifters
	firstdrifter <- 1
	lastdrifter <- 10000

	# The time record at which to make the plots
	time <- 6

	# The minimum and maximum longitudes and latitudes for the plots
	lonmin <- 11.6
	lonmax <- 27.4
	latmin <- -38.8
	latmax <- -27.8

	# The parameters below are used to plot the drifters density map
	# Number of grid cells in longitude direction
	nblon <- 100
	# Number of grid cells in latitude direction
	nblat <- 100
	# Number of grey levels
	nbgreylevels <- 50

	findnotdead <- function(x) seq(along=x)[x==0]

	nc <- open.ncdf(ncfile)

	# Checks for drifters that are 'not dead' i.e. still inside the domain
	nbdrifter <- lastdrifter-firstdrifter+1
	death <- get.var.ncdf(nc,'death',c(firstdrifter,time),c(nbdrifter,1))
	lon <- get.var.ncdf(nc,'lon',c(firstdrifter,time),c(nbdrifter,1))
	lat <- get.var.ncdf(nc,'lat',c(firstdrifter,time),c(nbdrifter,1))
	notdead <- findnotdead(death)
	lonnotdead <- lon[notdead]
	latnotdead <- lat[notdead]

	par(mfrow=c(2,2))

	# Plots locations of drifters as pixels
	lonlim <- c(lonmin,lonmax)
	latlim <- c(latmin,latmax)
	plot(lonnotdead,latnotdead,type='p',pch='.',xlim=lonlim,ylim=latlim,xlab='Longitude (°)',ylab='Latitude (°)')

	# Computes the number of drifters in each grid cell
	lontoplot <- seq(lonmin,lonmax,length=nblon)
	lattoplot <- seq(latmin,latmax,length=nblat)
	driftercount <- matrix(0,nrow=nblon,ncol=nblat)
	for (i in seq(1,length(notdead))){
		loni <- lonnotdead[i]
		lati <- latnotdead[i] 
		if ((loni>lonmin) & (loni<lonmax) & (lati>latmin) & (lati<latmax)){
			x <- 1+round((nblon-1)*(loni-lonmin)/(lonmax-lonmin))
			y <- 1+round((nblat-1)*(lati-latmin)/(latmax-latmin))
			driftercount[x,y] <- driftercount[x,y]+1
		}
	}

	colorbar <- gray((nbgreylevels:0)/nbgreylevels)

	# Plots drifters density map
	image(lontoplot,lattoplot,driftercount,col=colorbar,xlab='Longitude (°)',ylab='Latitude (°)')

	# Plots drifters density in prespective
	persp(lontoplot,lattoplot,driftercount,xlab='Longitude (°)',ylab='Latitude (°)',zlab='Number of drifters',theta=300,phi=30)

	# Plots the colorbar used in the drifters density map
	image(seq(0,max(driftercount),length=nbgreylevels),c(1),as.matrix(seq(1,nbgreylevels)),col=colorbar,xlab='Number of drifters',ylab='',yaxt='n')

	close.ncdf(nc)
}
