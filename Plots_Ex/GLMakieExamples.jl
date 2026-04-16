# GLMakie Examples
# https://github.com/JuliaPlots/Makie.jl
# https://makie.juliaplots.org/stable/

# %% Basic setup
using Random
using GLMakie
using Colors
using LinearAlgebra

# %% 1. Basic 2D Scatter Plot
println("1. Basic 2D Scatter Plot")
x = randn(100)
y = randn(100)

fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Basic 2D Scatter",
    xlabel = "x",
    ylabel = "y")

scatter!(ax, x, y, color = :blue, markersize = 8)
fig

# %% 2. Interactive 3D Scatter Plot
println("2. Interactive 3D Scatter Plot")
x_3d = randn(50)
y_3d = randn(50)
z_3d = randn(50)

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Interactive 3D Scatter",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

scatter!(ax, x_3d, y_3d, z_3d, 
    color = 1:50, 
    colormap = :viridis, 
    markersize = 15)
fig

# %% 3. 3D Surface Plot
println("3. 3D Surface Plot")
x_surf = -3:0.1:3
y_surf = -3:0.1:3
z_surf = [exp(-(x^2 + y^2)) for x in x_surf, y in y_surf]

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "3D Surface Plot",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

surface!(ax, x_surf, y_surf, z_surf, colormap = :plasma)
fig

# %% 4. 3D Line Plot
println("4. 3D Line Plot")
t = 0:0.01:4π
x_line = cos.(t)
y_line = sin.(t)
z_line = t ./ (4π)

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "3D Spiral",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

lines!(ax, x_line, y_line, z_line, 
    color = 1:length(t), 
    colormap = :rainbow, 
    linewidth = 3)
fig

# %% 5. Multiple 3D Objects
println("5. Multiple 3D Objects")
fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Multiple 3D Objects",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

# Sphere
sphere_points = 50
θ = range(0, 2π, sphere_points)
φ = range(0, π, sphere_points)
x_sphere = [cos(θi) * sin(φj) for θi in θ, φj in φ]
y_sphere = [sin(θi) * sin(φj) for θi in θ, φj in φ]
z_sphere = [cos(φj) for θi in θ, φj in φ]

surface!(ax, x_sphere, y_sphere, z_sphere, 
    color = :lightblue, 
    alpha = 0.7)

# Points
scatter!(ax, randn(20), randn(20), randn(20), 
    color = :red, 
    markersize = 10)

fig

# %% 6. Real-time Animation
println("6. Real-time Animation")
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Real-time Animation",
    xlabel = "x",
    ylabel = "y",
    limits = (-2, 2, -2, 2))

# Create a scatter plot that we'll update
scatter_obj = scatter!(ax, [0.0], [0.0], color = :red, markersize = 20)

# Animation function
function animate(frame)
    t = frame * 0.1
    x_pos = cos(t)
    y_pos = sin(t)
    
    # Update scatter plot data
    scatter_obj[1] = [x_pos]
    scatter_obj[2] = [y_pos]
    
    # Force redraw
    notify(scatter_obj)
end

# Start animation (this will run in the background)
# You can stop it by closing the window or pressing Ctrl+C
display(fig)

# %% 7. Interactive Volume Plot
println("7. Interactive Volume Plot")
# Create a 3D volume
vol_size = 32
volume_data = zeros(vol_size, vol_size, vol_size)

# Create a sphere in the volume
center = vol_size ÷ 2
radius = vol_size ÷ 4
for i in 1:vol_size, j in 1:vol_size, k in 1:vol_size
    dist = sqrt((i - center)^2 + (j - center)^2 + (k - center)^2)
    if dist < radius
        volume_data[i, j, k] = 1.0 - dist / radius
    end
end

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Interactive Volume Plot",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

volume!(ax, volume_data, 
    colormap = :viridis, 
    algorithm = :absorption)
fig

# %% 8. Mesh Plot
println("8. Mesh Plot")
using GeometryBasics

# Create a simple mesh
vertices = [
    Point3f(0, 0, 0),
    Point3f(1, 0, 0),
    Point3f(0, 1, 0),
    Point3f(1, 1, 0),
    Point3f(0, 0, 1),
    Point3f(1, 0, 1),
    Point3f(0, 1, 1),
    Point3f(1, 1, 1)
]

faces = [
    TriangleFace(1, 2, 3),
    TriangleFace(2, 4, 3),
    TriangleFace(5, 6, 7),
    TriangleFace(6, 8, 7),
    TriangleFace(1, 2, 5),
    TriangleFace(2, 6, 5),
    TriangleFace(3, 4, 7),
    TriangleFace(4, 8, 7),
    TriangleFace(1, 3, 5),
    TriangleFace(3, 7, 5),
    TriangleFace(2, 4, 6),
    TriangleFace(4, 8, 6)
]

mesh_obj = GeometryBasics.Mesh(vertices, faces)

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "Mesh Plot",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

mesh!(ax, mesh_obj, color = :orange)
fig

# %% 9. Contour Plot
println("9. Contour Plot")
x_contour = -3:0.1:3
y_contour = -3:0.1:3
z_contour = [sin(x) * cos(y) for x in x_contour, y in y_contour]

fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Contour Plot",
    xlabel = "x",
    ylabel = "y")

contour!(ax, x_contour, y_contour, z_contour, 
    colormap = :RdBu, 
    levels = 10)
fig

# %% 10. Interactive Scatter with Hover
println("10. Interactive Scatter with Hover")
x_hover = randn(20)
y_hover = randn(20)
labels = ["Point $i" for i in 1:20]

fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Interactive Scatter with Hover",
    xlabel = "x",
    ylabel = "y")

scatter!(ax, x_hover, y_hover, 
    color = 1:20, 
    colormap = :viridis, 
    markersize = 15)

# Add text labels that appear on hover (this is a simplified version)
for (i, (x, y, label)) in enumerate(zip(x_hover, y_hover, labels))
    text!(ax, x, y + 0.2, text = label, 
        align = (:center, :bottom), 
        textsize = 12)
end

fig

# %% 11. 3D Vector Field
println("11. 3D Vector Field")
x_vec = -2:0.5:2
y_vec = -2:0.5:2
z_vec = -2:0.5:2

# Create vector field data
u = zeros(length(x_vec), length(y_vec), length(z_vec))
v = zeros(length(x_vec), length(y_vec), length(z_vec))
w = zeros(length(x_vec), length(y_vec), length(z_vec))

for i in 1:length(x_vec), j in 1:length(y_vec), k in 1:length(z_vec)
    x, y, z = x_vec[i], y_vec[j], z_vec[k]
    r = sqrt(x^2 + y^2 + z^2)
    if r > 0.1
        u[i, j, k] = -y / r^2
        v[i, j, k] = x / r^2
        w[i, j, k] = 0.1
    end
end

fig = Figure()
ax = Axis3(fig[1, 1], 
    title = "3D Vector Field",
    xlabel = "x",
    ylabel = "y",
    zlabel = "z")

arrows!(ax, x_vec, y_vec, z_vec, u, v, w, 
    lengthscale = 0.3, 
    arrowcolor = :red, 
    arrowsize = 0.1)
fig

# %% 12. Real-time Data Visualization
println("12. Real-time Data Visualization")
fig = Figure()
ax = Axis(fig[1, 1], 
    title = "Real-time Data",
    xlabel = "Time",
    ylabel = "Value",
    limits = (0, 100, -2, 2))

# Initialize empty line plot
line_obj = lines!(ax, Float64[], Float64[], color = :blue, linewidth = 2)

# Data storage
time_data = Float64[]
value_data = Float64[]

# Update function for real-time data
function update_data()
    for i in 1:100
        push!(time_data, i)
        push!(value_data, sin(i * 0.1) + 0.1 * randn())
        
        # Update the line plot
        line_obj[1] = time_data
        line_obj[2] = value_data
        
        # Keep only last 100 points
        if length(time_data) > 100
            popfirst!(time_data)
            popfirst!(value_data)
        end
        
        notify(line_obj)
        sleep(0.1)  # Update every 100ms
    end
end

# Start the real-time update
display(fig)

println("All GLMakie examples completed!")
println("Note: Some interactive features require the plot window to be open.") 