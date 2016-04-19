
library(ncdf4)

setwd("~/Dropbox/GitHub/pmca/prohaska")

nc <- nc_open("slp.mnmean.nc")   # opens nc file
slp.lon <- ncvar_get( nc, "lon")
slp.lat <- ncvar_get( nc, "lat")
slp.t <- ncvar_get( nc, "time")
slp.raw <- ncvar_get(nc, "slp")
nc_close(nc)
