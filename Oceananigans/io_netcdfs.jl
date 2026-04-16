##
using NetCDF
using Plots

file = "gom_t009.nc"
println(ncinfo(file))

##
data = ncread(file, "water_temp")
println(size(data))

heatmap(data[:,:,1,1], colorbar=false, title="Temp")])
