# sudo apt-get install gmt gmt-dcw gmt-gshhg ghostscript
# optional 
# sudo apt-get install gdal-bin graphicsmagick ffmpeg 
# Make link sudo ln -s /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so
# READ!!!!!! If error with GLIBC then do cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/julia/lib/julia/
using GMT
using NCDatasets

## ---------------- READING DATA --------------------------------
file_name = "/data/UN_Litter_data/output/YesWinds_YesDiffusion_NoUnbeaching_2010_01.nc"
ds = Dataset(file_name)
# Assuming we don't have missing valuesj
lat =  ds["lat"]
lon =  ds["lon"]
clon = 180

timesteps = size(lat)[1]
# Threads.@threads for i in 1:10
for i in 1:Int(365*2)
# for i in 1:2
    # println("$i --> thread id $(Threads.threadid())")

    # gmthelp(coast) # if you want to see what each option means
    # Simple white background
    # res=full, high, intermediate, low or crude
    # coast(R=[0+clon 360+clon -60 76], J=:Mercator, res=:low, land=:grey, borders=(a=:1), C=:blue, title="World litter day $i")
    # With bathymetri as a background
    grdimage("@earth_relief_20m.grd",  colorbar=true, R=[0+clon 360+clon -60 76], J=:Mercator, 
    coast=( res=:crude, borders=(a=:1), C=:blue, title="World litter day $i"))

    scatter!(convert(Array{Float32},lon[i,:]), convert(Array{Float32},lat[i,:]), fmt=:png, marker=:circle, 
                width=1, size=0.03, markerfacecolor=:red, name="imgs/$(lpad(i,4,'0')).png")
end


##
grdcontour("@HI_topo_02.nc", cont=1, show=true)