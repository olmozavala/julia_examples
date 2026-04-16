##
using Plots
using Oceananigans
using NetCDF

##
test_fu(x, y, z) = .5
test_fv(x, y, z) = 0
u_forcing = Forcing(test_fu)
v_forcing = Forcing(test_fv)

grid = RectilinearGrid(CPU(), size=(128, 128), x=(0, 2π), y=(0, 2π), topology=(Periodic, Periodic, Flat))
model = NonhydrostaticModel(; grid, forcing=(u=u_forcing, v=v_forcing), advection=WENO())  # Working
# model = HydrostaticFreeSurfaceModel(; grid, forcing=(u=u_forcing, v=v_forcing), advection=WENO()) # Not working
# model = Hy(; grid, forcing=(u=u_forcing, v=v_forcing), advection=WENO())
                                                        
# ϵ(x, y) = 0
ϵ(x, y) = 2rand() - 1

set!(model, u=ϵ, v=ϵ)
simulation = Simulation(model; Δt=0.01, stop_time=4)

##
output_file = "helloworld.nc"
fields = Dict("u" => model.velocities.u, "v" => model.velocities.v)

simulation.output_writers[:growth] = NetCDFOutputWriter(model, fields, 
                                                        filename = output_file,
                                                        schedule = TimeInterval(.1),
                                                        overwrite_existing = true)

##
run!(simulation)
println("Done!")
##
println(ncinfo(output_file))
data = ncread(output_file, "u")
println(size(data))
times = size(data)[4]
##
anim = @animate for i = 1:times
    # Your heatmap data generation or retrieval logic here
    # For example, let's say you have a function `get_heatmap_data(i)`
    # which returns the data for the i-th heatmap
    heatmap(data[:,:,1,i], colorbar=false, title="u time $i")
end
gif(anim, "./heatmap_animation.gif", fps = 5)

##
println(size(model.velocities.v))
# println(model.velocities.u[:,:,1])

##
v = interior(model.velocities.v,:,:,1)
u = interior(model.velocities.u,:,:,1)
z = znodes(model.velocities.v)
println(typeof(z))

##
l = @layout[a b]
p1 = heatmap(Array(u))
p2 = heatmap(Array(interior(model.velocities.v,:,:,1)))
plot(p1, p2, layout=l)

## Plotting mulitple times
# @gif for i=1:1500
#     step!(attractor)
#     push!(plt, attractor.x, attractor.y, attractor.z)
# end every 10
