plot_traj_Ichthyop_v3_ggmap <- function(){

	library(ncdf)
	library(ggmap)

	# The netcdf input file
        ncfile <- '/home/clett/IRD/Ichthyop/Ichthyop3.2b/ichthyop-3.2b/ichthyop-3.2b/output/roms3d_ichthyop-run201509081630.nc'
	
	# In case one wishes to consider only a subset of all drifters
	# Index of the first and last drifters
        firstdrifter <- 1
	lastdrifter <- 40
	
	# In case one wishes to consider only a subset of all time records
	# Index of the first and last time records
	firsttime <- 1
	lasttime <- 61

	# The minimum and maximum longitudes and latitudes for the plots    
        lonmin <- 15.0
	lonmax <- 20.0
	latmin <- -35.0
	latmax <- -30.0
	
	# Number of colors
	nbcolor <- 50
	
	nc <- open.ncdf(ncfile)
        nbdrifter <- lastdrifter - firstdrifter + 1
        nbtime <- lasttime - firsttime + 1
	drifter <- rep(seq(firstdrifter, lastdrifter), each = nbtime)
        lon <- as.vector(t(get.var.ncdf(nc, 'lon', c(firstdrifter, firsttime), c(nbdrifter, nbtime))))
        lat <- as.vector(t(get.var.ncdf(nc, 'lat', c(firstdrifter, firsttime), c(nbdrifter, nbtime))))
        depth <- as.vector(t(get.var.ncdf(nc, 'depth', c(firstdrifter, firsttime), c(nbdrifter, nbtime))))
        df <- data.frame(drifter, lon, lat, depth)
        close.ncdf(nc)
        
        lonlim <- c(lonmin, lonmax)
	latlim <- c(latmin, latmax)
        mymap <- get_map(location = c(lon = (lonmin + lonmax) / 2, lat = (latmin + latmax) / 2), zoom = 7, maptype = "hybrid")
	map <- ggmap(mymap)
        map <- map + geom_point(data = df, aes(x = lon, y = lat, colour = drifter), size = 1) + scale_colour_gradientn(colours = rainbow(nbcolor))
        map + labs(x = 'Longitude (°E)', y = 'Latitude (°N)')
}