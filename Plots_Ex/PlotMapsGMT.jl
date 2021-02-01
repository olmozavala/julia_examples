# Install on Ubuntu
# sudo apt-get install gmt gmt-dcw gmt-gshhg ghostscript
# optional 
# sudo apt-get install gdal-bin graphicsmagick ffmpeg 
# Make link sudo ln -s /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so
# READ!!!!!! If error with GLIBC then do cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/julia/lib/julia/
using GMT

# Examples at: https://www.generic-mapping-tools.org/GMT.jl/latest/examples/

## Hello World
plot(1:10, rand(10), lw=1, lc=:blue, fmt=:png, marker=:square,
     markeredgecolor=0, size=0.1, markerfacecolor=:red, title="Hello World",
     xlabel="Spoons", ylabel="Forks", show=true)

## Hello earth
x = range(0, stop=2pi, length=180);       seno = sin.(x/0.2)*45;
coast(region=[0 360 -90 90], proj=(name=:laea, center=(300,30)), frame=:g,
      res=:crude, land=:navy, figsize=6)

plot!(collect(x)*60, seno, lw=0.5, lc=:red, fmt=:png, marker=:circle,
      markeredgecolor=0, size=0.05, markerfacecolor=:cyan, show=true)


## Color images
topo = makecpt(color=:rainbow, range=(1000,5000,500), continuous=true); # Creates a color palette
grdimage("@tut_relief.nc", shade=(azimuth=100, norm="e0.8"), proj=:Mercator, frame=:a, color=topo)
# Assign a colorbar
colorbar!(pos=(anchor=:TC,length=(12.5,0.6), horizontal=true, offset=(0,1.0)),
          color=topo, frame=(ylabel=:m,), fmt=:jpg, show=true)


## Persperctive view
topo = makecpt(color=:rainbow, range=(1000,5000,500), continuous=true);
grdview("@tut_relief.nc", proj=:Mercator, zsize=1, shade=(azim=100, norm="e0.8"), view=(135,30),
        frame=:a, fmt=:jpg, Q="i100", show=true)

##  Contours
G = GMT.peaks();
grdcontour(G, cont=1, annot=2, fmt=:png, show=true, scale=.5)

## My own sine contour
using LinearAlgebra
x = sin.(1:.1:2*π)
s = size(x)[1]
M = ones(s,s)
A = M.*x + transpose(M.*x) 
grdcontour(A, cont=1, annot=2, fmt=:png, show=true, scale=.5)


## Image from gridded data
using NCDatasets
file_name = "/home/olmozavala/Dropbox/TestData/GIS/Bathymetry/Bathymetry.nc"
ds = Dataset(file_name)
data = permutedims(ds["Bathymetry"][:])
grdimage(data, region="0,360,-90,90", show=true, portrait=false)

## My template
coast(region=[-176 -4 -15 50],                                   # The Map limits	
      proj=(name=:Albers, center=[-83 22], parallels=[25 45]),  # The projection parameters
      frame=:ag,          # Tell it to set annotations and grid lines automatically
      resolution=:low,    # Use the low resolution coastlines
      area=250,           # Do not plot polygons with areas < 250 km^2
      land=:green,        # Paint land with green
      shore=:thinnest,    # Coastlines are drwan with a 0.1 pt thickness
      fmt=:png,
      show=true)
    
