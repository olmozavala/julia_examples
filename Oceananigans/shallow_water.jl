##
using Oceananigans
using Oceananigans.Models: ShallowWaterModel

grid = RectilinearGrid(GPU(), size = (48, 128),
                       x = (0, 2π),
                       y = (-10, 10),
                       topology = (Periodic, Bounded, Flat))

gravitational_acceleration = 1
coriolis = FPlane(f=1)
##

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                            timestepper = :RungeKutta3,
                            momentum_advection = WENO())


##
U = 1  # Maximum jet velocity
H = 10 # Reference depth
f = coriolis.f
g = gravitational_acceleration
Δη = f * U / g  # Maximum free-surface deformation as dictated by geostrophy

h̄(x, y) = H - Δη * tanh(y)
ū(x, y) = U * sech(y)^2

ω̄(x, y) = 2 * U * sech(y)^2 * tanh(y)

small_amplitude = 1e-4

 uⁱ(x, y) = ū(x, y) + small_amplitude * exp(-y^2) * randn()
uhⁱ(x, y) = uⁱ(x, y) * h̄(x, y)

ū̄h(x, y) = ū(x, y) * h̄(x, y)

set!(model, uh = ū̄h, h = h̄)
##

uh, vh, h = model.solution

# Build velocities
u = uh / h
v = vh / h

# Build and compute mean vorticity discretely
ω = Field(∂x(v) - ∂y(u))
compute!(ω)

# Copy mean vorticity to a new field
ωⁱ = Field((Face, Face, Nothing), model.grid)
ωⁱ .= ω

# Use this new field to compute the perturbation vorticity
ω′ = Field(ω - ωⁱ)

##
simulation = Simulation(model, Δt = 1e-2, stop_time = 100)

## 
using LinearAlgebra: norm

perturbation_norm(args...) = norm(v)

growth_filename = joinpath(@__DIR__, "shallow_water_Bickley_jet_perturbation_norm.nc")
simulation.output_writers[:growth] = NetCDFOutputWriter(model, (; perturbation_norm),
                                                        filename = growth_filename,
                                                        schedule = IterationInterval(1),
                                                        dimensions = (; perturbation_norm = ()),
                                                        overwrite_existing = true)

run!(simulation)

