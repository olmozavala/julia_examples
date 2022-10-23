using Missings

"""
This function replaces an Array{Union{Missing, Any}} with Array{Float64}
"""
function missingsToFloat(data, repmissing = 0.0)
    return collect(Missings.replace(data,repmissing))
end

"""
This function obtains the bbox for a Dataset. It assumes the names for the 
    coordinates are LAT and LON
"""
function getBBOX(ds)
    field_names = keys(ds)

    if ("LAT" in field_names) & ("LON" in field_names) 
        lats = ds["LAT"][:]
        lons = ds["LON"][:]
        bbox = (minimum(lons), maximum(lons), minimum(lats), maximum(lats))
        return bbox
    else
        return -1
    end
end