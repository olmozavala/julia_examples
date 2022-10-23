# Install on Ubuntu
# sudo apt-get install gmt gmt-dcw gmt-gshhg ghostscript
# optional 
# sudo apt-get install gdal-bin graphicsmagick ffmpeg 
# Make link sudo ln -s /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so
# READ!!!!!! If error with GLIBC then do cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/julia/lib/julia/
using GMT

# General info 
# Modules summary: https://www.generic-mapping-tools.org/GMT.jl/dev/modules/#By-Modules
# Common options: https://www.generic-mapping-tools.org/GMT.jl/dev/common_opts/
# Projections: https://docs.generic-mapping-tools.org/latest/cookbook/map-projections.html?highlight=projections
# Examples at: https://www.generic-mapping-tools.org/GMT.jl/latest/examples/

## --------------- Hello World -----------
plot(1:10, rand(10), lw=1, lc=:blue, fmt=:png, marker=:square,
     markeredgecolor=0, size=0.1, markerfacecolor=:red, title="Hello World",
     xlabel="Spoons", ylabel="Forks", show=true)

## -------------  Hello earth ------------
x = range(0, stop=2pi, length=180)
seno = sin.(x/0.2)*45

coast(region=[0 360 -30 30], proj=:merctor, frame=:g,
      res=:crude, land=:navy, figsize=6, show=true, fmt=:png)

plot!(collect(x)*60, seno, lw=0.5, lc=:red, fmt=:png, marker=:circle,
      markeredgecolor=0, size=0.05, markerfacecolor=:cyan, show=true)

## Subplots
subplot(grid=(2,2), panels_size=8, region=(0, 100, 0, 80), margins="5p", autolabel=true, col_axes=(bott=true,), 
row_axes=(left=true,), axes="wstr", name="panels.pdf")
subplot(:set)
      basemap(region=(0,80,0,50))
      subplot(:set)
      plot([50 40], marker=:circle, mc=:red)
      subplot(:set)
      plot([50 40], marker=:square, mc=:green)
      subplot(:set, panel=(2,2))
      plot([50 40], marker=:star, mc=:blue)
subplot(show=true, fmt=:png)

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
grdimage(data,  show=true, portrait=false,  M=true, fmt=:png)
grdimage("@earth_relief_20m.grd",  colorbar=true, coast=true, show=true)


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
    
## Arrows
# x,y, w, l
# arrows([0 8.2 0 6], limits=(-1,4,7,9), arrow=(len=2,start=:arrow,stop=:tail,shape=0.5),

n = 10
l = 4
A = [rand(n).*10 rand(n).*10 rand(n).*360 ones(n).*l] # x, y, angle, length
# arrows([0 8.2 90 6], limits=(minimum(A[:,1]), maximum(A[:,1]), -10,20), arrow=(len=2,start=:arrow,stop=:tail,shape=0.1),
#[x, y, angle, length]  and limits (xmin, xmax, ymin, ymax)
bbox = (minimum(A[:,1])-l, maximum(A[:,1])+l, minimum(A[:,2])-l, maximum(A[:,2])+l)
arrows(A, limits=bbox, arrow=(len=1,start=:tail,stop=:arrow,shape=.2),
            figsize=(12,4), axis=(grid=1.0, annot=1), pen="2p", show=true)
##
# for x = 0:5
#       arrows([x 8.2 0 3], limits=(-1,4,7,9), arrow=(len=1,start=:arrow,stop=:tail,shape=0.5),
#            proj=:X, figsize=(12,4), axis=:a, pen="4p", show=true)
#       sleep(1)
# end

# arrows([0.0 0 5 5], limits=(0,5,0,5), proj=:X10, axis=(annot=0.5, ticks=0.25, grid=0.5),
#           arrow=(len=0.5,stop=1,uv=0.5), show=true)


# arrows([0.0 0 5 5], limits=(0,5,0,5), proj=:X10, axis=(annot=0.5, ticks=0.25, grid=0.5),
#           arrow=(len=0.5,stop=1,uv=0.5), show=true)
