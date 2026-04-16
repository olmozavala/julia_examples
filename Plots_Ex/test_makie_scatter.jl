# Test script to verify Makie scatter functions work correctly
using CairoMakie

println("Testing Makie scatter functions...")

# Test 1: Basic scatter
x = 1:10
y = rand(10)

fig = Figure()
ax = Axis(fig[1, 1], title = "Test Scatter")
scatter!(ax, x, y, color = :blue, markersize = 8)
fig

println("✓ Basic scatter test passed")

# Test 2: Scatter with color values (using colormap)
colors = rand(10)  # Single values for colormap
fig2 = Figure()
ax2 = Axis(fig2[1, 1], title = "Test Scatter with Color Values")
scatter!(ax2, x, y, color = colors, markersize = 8, colormap = :viridis)
fig2

println("✓ Scatter with color values test passed")

# Test 3: Scatter with RGB colors (proper way)
using Colors
rgb_colors = [RGB(rand(), rand(), rand()) for _ in 1:10]
fig3 = Figure()
ax3 = Axis(fig3[1, 1], title = "Test Scatter with RGB Colors")
scatter!(ax3, x, y, color = rgb_colors, markersize = 8)
fig3

println("✓ Scatter with RGB colors test passed")

# Test 4: Scatter with points
points = Point2f.(x, y .+ 1)
fig4 = Figure()
ax4 = Axis(fig4[1, 1], title = "Test Scatter with Points")
scatter!(ax4, points, color = 1:10, markersize = 8, colormap = :viridis)
fig4

println("✓ Scatter with points test passed")

println("All scatter tests completed successfully!") 