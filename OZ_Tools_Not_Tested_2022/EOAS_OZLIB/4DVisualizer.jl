# using Plots
using NCDatasets
using GMT
include("/home/olmozavala/Dropbox/TutorialsByMe/Julia/Julia_Examples/OZ_Tools/EOAS_OZLIB/GeoUtils.jl")

"""
It generates plots 
"""
# function plotSingleDepthFieldsFromDS(ds, field_names=[], file_prefix="", fig_per_col=1)
function plotSingleDepthFieldsFromDS(ds, field_names=[], fig_per_col=1)
    if isempty(field_names)
        field_names = keys(ds)
    end

    tot_fields = length(field_names)
    tot_rows = Int(ceil(tot_fields/fig_per_col))
    bbox = getBBOX(ds)

    # This part is working
    for field_name in field_names
        tic()
        t = missingsToFloat(ds[field_name]) # TODO line is too slow, if I don't make a separate expression it won't load
        toc()
        grdimage(permutedims(t))
        coast!(region=bbox,# The Map limits	
            proj=(name=:mercator),  # The projection parameters
            frame=:ag,          # Tell it to set annotations and grid lines automatically
            resolution=:coars,    # Use the low resolution coastlines
            area=250,           # Do not plot polygons with areas < 250 km^2
            shore=:thinnest,    # Coastlines are drwan with a 0.1 pt thickness
            land=:grey, savefig="$(field_name).png", borders=(a=:1), C=:blue, title="$field_name")
    end

    # subplot(grid=(tot_rows, fig_per_col), region=bbox, panels_size=8, margins="5p", col_axes=(bott=true,))
    # for field_name in field_names
    #     println("Plotting $field_name....")
    #     field = missingsToFloat(ds[field_name])

    #     subplot(:set)
    #     grdimage(permutedims(field))
    #     coast!(                                 # The Map limits	
    #         proj=(name=:mercator),  # The projection parameters
    #         frame=:ag,          # Tell it to set annotations and grid lines automatically
    #         resolution=:low,    # Use the low resolution coastlines
    #         area=250,           # Do not plot polygons with areas < 250 km^2
    #         shore=:thinnest,    # Coastlines are drwan with a 0.1 pt thickness
    #         land=:grey,
    #         borders=(a=:1),
    #         C=:blue,
    #         title="$field_name")
    # end
    # subplot(:show)

    println("Done!")
end