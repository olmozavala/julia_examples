using NCDatasets
using GMT
include("/home/olmozavala/Dropbox/TutorialsByMe/Julia/Julia_Examples/OZ_Tools/EOAS_OZLIB/4DVisualizer.jl")
include("/home/olmozavala/Dropbox/TutorialsByMe/Julia/Julia_Examples/OZ_Tools/EOAS_OZLIB/GeoUtils.jl")
##

# folder =  "/data/HYCOM/DA_HYCOM_TSIS/preproc"
folder =  "/home/olmozavala/Dropbox/TestData/netCDF/HYCOM/S_T_U_V_Surface"
all_files = readdir(folder)
files = [x for x in all_files if occursin("model", x) ]# Read folders 

# Plot variables from single file
ds = Dataset(joinpath(folder, files[1]))
all_fields = keys(ds)
bbox = getBBOX(ds)

lats = ds["LAT"][:]
lons = ds["LON"][:]
##
plotSingleDepthFieldsFromDS(ds, ["LAT","temp", "srfhgt"], 2)

##
subplot(grid=(2,2), panels_size=8, region=(0, 100, 0, 80), margins="5p", col_axes=(bott=true,), row_axes=(left=true,), axes="wstr", name="panels.pdf")
subplot(:set)
coast(region=bbox, proj=(name=:mercator), frame=:ag, resolution=:low, area=250, shore=:thinnest, land=:grey, borders=(a=:1), C=:blue, fmt=:png, show=true)
subplot(:set)
coast(region=bbox, proj=(name=:mercator), frame=:ag, resolution=:low, area=250, shore=:thinnest, land=:white, borders=(a=:1), C=:blue, fmt=:png, show=true)
subplot(:set)
coast(region=bbox, proj=(name=:mercator), frame=:ag, resolution=:low, area=250, shore=:thinnest, land=:red, borders=(a=:1), C=:blue, fmt=:png, show=true)
subplot(:set)
coast(region=bbox, proj=(name=:mercator), frame=:ag, resolution=:low, area=250, shore=:thinnest, land=:blue, borders=(a=:1), C=:blue, fmt=:png, show=true)
subplot(:show)