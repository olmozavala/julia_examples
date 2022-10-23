using GeoArrays
using Plots
# https://github.com/evetion/GeoArrays.jl

## Read and plot example file
root = "/home/olmozavala/Dropbox/TestData/GIS/Topo/"
bath = GeoArrays.read(joinpath(root,"merged.tif"))
println("Done reading!")
# Access data with bath.A

## Useful stuff
println(bath.crs) # Acess projection
# println(coords(bath,[1,1]))  # Accessing coordinated in specific location
println(indices(bath,[0.0,0.0]))  # Accessing indices in specific location
heatmap(bath, band=1)

## ----------------------- Save data to netcdf ---------------
using NCDatasets
# The mode "c" stands for creating a new file (clobber)
ds = Dataset(joinpath(root,"Bathymetry.nc"), "c")
dims = size(bath.A)

# Define the dimension "lon" and "lat" with the size 100 and 110 resp.
defDim(ds,"lon",dims[1])
defDim(ds,"lat",dims[2])

# Define a global attribute
ds.attrib["title"] = "Global Bathymetry"

# Define the variables temperature with the attribute units
v = defVar(ds,"Bathymetry",Float32,("lon","lat"), attrib = Dict( "units" => "pixels"))

# add additional attributes
v.attrib["comments"] = "this is a string attribute with Unicode Ω"

v[:,:] = bath.A[:,end:-1:1,1]
close(ds)

## SRTM dat https://topex.ucsd.edu/WWW_html/srtm15_plus.html
# https://www.ngdc.noaa.gov/mgg/topo/gltiles.html
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/