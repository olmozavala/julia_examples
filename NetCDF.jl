# ************************ USING NCDatasets (Prefered) *************************
using NCDatasets
using Makie
using Images  # ImageMagick ImageIO FileIO
using FileIO
using LinearAlgebra

## ---------------- READING DATA --------------------------------
file_name = "/home/olmozavala/Dropbox/TestData/netCDF/gfs.nc"
ds = Dataset(file_name)
print(ds)

# ----- Read and plotting  a field -----
## Single field
field = ds["TMP_P0_L104_GLL0"]
nfield = transpose(field./maximum(field))

# imgg = Gray.(nfield)

## Multiple fields
fields_to_plot = ["TMP_P0_L7_GLL0", "TMP_P0_L104_GLL0", "RH_P0_L204_GLL0"]

for c_field_name in fields_to_plot
    field = ds[c_field_name]
    println("*********** Drawing field: $c_field_name with size: $(size(field))*********** ")
    println("Attributes: $(field.attrib)")
    println("Type of field: $(typeof(field))")
    # push!(all_scenes, title(heatmap(reverse(field[:,:], dims=2)), c_field_name))
end
close(ds)
## ---------------- WRITING DATA --------------------------------

# The mode "c" stands for creating a new file (clobber)
ds = Dataset("test.nc","c")

# Define the dimension "lon" and "lat" with the size 100 and 110 resp.
defDim(ds,"lon",100)
defDim(ds,"lat",110)

# Define a global attribute
ds.attrib["title"] = "this is a test file"

# Define the variables temperature with the attribute units
v = defVar(ds,"temperature",Float32,("lon","lat"), attrib = OrderedDict( "units" => "degree Celsius"))

# add additional attributes
v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

# Generate some example data
data = [Float32(i+j) for i = 1:100, j = 1:110]

# write a single column
v[:,1] = data[:,1]

# write a the complete data set
v[:,:] = data

close(ds)

## ************************ USING NetCDF *************************
using NetCDF
using Colors
using Makie
using AbstractPlotting: hbox, vbox
using AbstractPlotting

# ---------------- READING DATA --------------------------------
# file_name = "/home/olmozavala/Dropbox/TestData/netCDF/HYCOM/hycom_GLBv0.08_536_2010010112_t000.nc"
file_name = "/home/olmozavala/Dropbox/TestData/netCDF/gfs.nc"

## ------ Summary of file -------
print(ncinfo(file_name))

## ----- Read and plotting  a field -----
fields_to_plot = ["TMP_P0_L7_GLL0", "TMP_P0_L104_GLL0", "RH_P0_L204_GLL0"]
scene = Scene() # Makies parent scene

all_scenes = Scene[]
for c_field_name in fields_to_plot
    println("Drawing field: $c_field_name")
    x  = NetCDF.open(file_name, c_field_name);
    println("Size of field: $(size(x))")
    push!(all_scenes, title(heatmap(reverse(x, dims=2)), c_field_name))
end

s = hbox(all_scenes)
## ----- Read coordinates

# ---------------- WRITING DATA --------------------------------

## ---------------- Example to plot bathymetry --------------------------------
file_name = "/home/olmozavala/Dropbox/TestData/GIS/Bathymetry/gebco_1min_-98_18_-78_31.nc"
ds = Dataset(file_name)
print(ds)