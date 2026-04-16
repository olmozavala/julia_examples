# Simple GLMakie Examples
# GLMakie is great for interactive 3D plots and real-time graphics

# %% Basic setup
using GLMakie
using Random

# %% 1. Simple 2D Plot
println("1. Simple 2D Plot")
x = 0:0.1:2π
y = sin.(x)

fig = Figure()
ax = Axis(fig[1, 1], title = "Simple Sine Wave")
lines!(ax, x, y, color = :blue, linewidth = 2)
fig

# %% 2. 3D Scatter Plot
println("2. 3D Scatter Plot")
x = randn(100)
y = randn(100)
z = randn(100)

fig = Figure()
ax = Axis3(fig[1, 1], title = "3D Scatter")
scatter!(ax, x, y, z, color = 1:100, markersize = 10)
fig

# %% 3. 3D Surface Plot
println("3. 3D Surface Plot")
x = -2:0.1:2
y = -2:0.1:2
z = [exp(-(xi^2 + yi^2)) for xi in x, yi in y]

fig = Figure()
ax = Axis3(fig[1, 1], title = "3D Surface")
surface!(ax, x, y, z, colormap = :viridis)
fig

# %% 4. 3D Line Plot
println("4. 3D Line Plot")
t = 0:0.01:4π
x = cos.(t)
y = sin.(t)
z = t ./ (4π)

fig = Figure()
ax = Axis3(fig[1, 1], title = "3D Spiral")
lines!(ax, x, y, z, color = 1:length(t), linewidth = 2)
fig

# %% 5. Multiple Plots in One Figure
println("5. Multiple Plots in One Figure")
fig = Figure()

# 2D plot
ax1 = Axis(fig[1, 1], title = "2D Plot")
x = 0:0.1:2π
lines!(ax1, x, sin.(x), color = :blue, label = "sin(x)")
lines!(ax1, x, cos.(x), color = :red, label = "cos(x)")
axislegend(ax1)

# 3D plot
ax2 = Axis3(fig[1, 2], title = "3D Plot")
scatter!(ax2, randn(50), randn(50), randn(50), color = :green, markersize = 8)

fig

# %% 6. Interactive Scatter with Colors
println("6. Interactive Scatter with Colors")
x = randn(50)
y = randn(50)
colors = rand(50)  # Values for colormap

fig = Figure()
ax = Axis(fig[1, 1], title = "Colored Scatter")
scatter!(ax, x, y, color = colors, markersize = 15, colormap = :plasma)
fig

# %% 7. Simple Animation
println("7. Simple Animation")
fig = Figure()
ax = Axis(fig[1, 1], title = "Animated Point", limits = (-2, 2, -2, 2))

# Create a scatter plot
scatter_obj = scatter!(ax, [0.0], [0.0], color = :red, markersize = 20)

# Animation function
function animate(frame)
    t = frame * 0.1
    x_pos = cos(t)
    y_pos = sin(t)
    
    # Update position
    scatter_obj[1] = [x_pos]
    scatter_obj[2] = [y_pos]
    
    # Force redraw
    notify(scatter_obj)
end

# Display the figure (animation will run)
display(fig)

println("All simple GLMakie examples completed!")
println("Try rotating and zooming the 3D plots with your mouse!") 