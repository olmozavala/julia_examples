# Julia Examples Collection

A comprehensive collection of self-contained Julia examples covering various domains including scientific computing, machine learning, differential equations, plotting, and more.

## 🖥️ Installing Julia on Linux

### **Recommended: Juliaup Installation**

Juliaup is the official Julia installer that manages multiple Julia versions efficiently.

#### **Install Juliaup:**
```bash
curl -fsSL https://install.julialang.org | sh
```

#### **Verify Installation:**
```bash
julia --version
```

## 📚 Package Management Best Practices

### **1. Project-Based Workflow (Recommended)**
Each example folder uses Julia's project system for dependency isolation:

```julia
# Navigate to project directory
cd("path/to/example")

# Start Julia REPL
julia

# Activate the project environment
using Pkg
Pkg.activate(".")

# Install dependencies (if needed)
Pkg.instantiate()

# Run the example
include("example.jl")
```

### **2. Adding Packages to a Project**
```julia
using Pkg

# Add a single package
Pkg.add("PackageName")

# Add multiple packages
Pkg.add(["Package1", "Package2", "Package3"])

# Add with specific version
Pkg.add(name="PackageName", version="1.2.3")

# Add from GitHub
Pkg.add(url="https://github.com/user/repo")
```

### **3. Package Development Workflow**
```julia
# Enter package mode
using Pkg
Pkg.activate(".")

# Add packages for development
Pkg.add("PackageName")

# Update all packages
Pkg.update()

# Remove a package
Pkg.rm("PackageName")

# Show project status
Pkg.status()

# Show package dependencies
Pkg.dependencies()
```

### **4. Environment Management**
```julia
# Create a new project
Pkg.generate("MyNewProject")

# Activate existing project
Pkg.activate("path/to/project")

# Activate global environment (not recommended for examples)
Pkg.activate()

# Show active environment
Pkg.project()
```

### **5. Package Registry Management**
```julia
# List available registries
Pkg.Registry.status()

# Add a registry
Pkg.Registry.add("RegistryName")

# Remove a registry
Pkg.Registry.rm("RegistryName")
```

## 🔧 Development Environment Setup

### **VS Code with Julia Extension (Recommended)**
1. Install [VS Code](https://code.visualstudio.com/)
2. Install Julia extension: `julialang.language-julia`
3. Configure Julia path in settings
4. Features: IntelliSense, debugging, plotting, REPL integration

### **Jupyter Notebooks**
```julia
# Install IJulia
using Pkg
Pkg.add("IJulia")

# Launch Jupyter
using IJulia
notebook()
```

### **Pluto.jl (Reactive Notebooks)**
```julia
# Install Pluto
using Pkg
Pkg.add("Pluto")

# Launch Pluto
using Pluto
Pluto.run()
```

## ⚡ Performance Tips

### **1. Package Precompilation**
```julia
# Force precompilation of all packages
using Pkg
Pkg.precompile()

# Precompile specific package
using PackageName  # This triggers precompilation
```

### **2. Startup Optimization**
Create `~/.julia/config/startup.jl`:
```julia
# Enable threading
ENV["JULIA_NUM_THREADS"] = "auto"

# Set number of threads manually
ENV["JULIA_NUM_THREADS"] = "4"

# Optimize for performance
ENV["JULIA_OPTIMIZE"] = "3"
```

### **3. Package Loading Optimization**
```julia
# Use Revise for development
using Pkg
Pkg.add("Revise")
using Revise

# Lazy loading with Requires.jl
using Pkg
Pkg.add("Requires")
```

## 🐛 Troubleshooting

### **Common Issues:**

1. **Package not found:**
   ```julia
   Pkg.update()  # Update package registry
   Pkg.add("PackageName")  # Try adding again
   ```

2. **Version conflicts:**
   ```julia
   Pkg.resolve()  # Resolve dependency conflicts
   Pkg.status()   # Check current status
   ```

3. **Precompilation errors:**
   ```julia
   Pkg.precompile()  # Force recompilation
   # Or delete ~/.julia/compiled and restart Julia
   ```

4. **Memory issues:**
   ```julia
   # Increase memory limit
   ENV["JULIA_MEMORY_LIMIT"] = "4G"
   ```

### **Getting Help:**
- [Julia Documentation](https://docs.julialang.org/)
- [Pkg.jl Documentation](https://pkgdocs.julialang.org/)
- [Julia Discourse](https://discourse.julialang.org/)
- [Stack Overflow](https://stackoverflow.com/questions/tagged/julia)

## 🚀 Getting Started

Each folder contains self-contained examples with their own `Project.toml` and `Manifest.toml` files. To run examples:

1. Navigate to the desired folder
2. Start Julia in that directory
3. Activate the project: `using Pkg; Pkg.activate(".")`
4. Run the example files

## 📁 Directory Structure

### 🎯 **Basics/**
Core Julia programming fundamentals and syntax examples.
- `Basics.jl` - Basic Julia syntax and operations
- `Arrays_oz.jl` - Array manipulation and operations
- `Strings_oz.jl` - String handling and text processing
- `Functions_oz.jl` - Function definition and usage
- `Loops_Ifs_Whiles_Compren_etc.jl` - Control structures and comprehensions
- `MyModule.jl` - Module creation and organization
- `MyFunctionsFile.jl` - Custom function definitions

### 🧠 **FLUX_Ex/**
Neural network and machine learning examples using Flux.jl.
- `0_Zygote.jl` - Automatic differentiation with Zygote
- `1_Flux_1_Hello_World.jl` - Basic Flux neural network example
- `2_Custom_Layer.jl` - Creating custom neural network layers
- `2_Custom_Loss.jl` - Implementing custom loss functions
- `OLD/` - Legacy examples and experiments
- `NotWorking/` - Experimental code that needs debugging

### 🔬 **DifferentialEquationsjl/**
Ordinary and partial differential equation solving examples.
- `1_ODE_Hello_World.jl` - Basic ODE solving introduction
- `2_ODE_Lorenz96.jl` - Lorenz 96 system implementation
- `3_ODE_Examples_Real_life.jl` - Real-world ODE applications
- `NotWorking/` - Experimental differential equation code

### 🤖 **SciML_Ex/**
Scientific Machine Learning (SciML) examples and applications.
- **NODEs/** - Neural Ordinary Differential Equations
  - `0_ODE_Hello_World.jl` - ODE basics for NODEs
  - `1_NeuralODE.jl` - Neural ODE implementation
- **SINDY/** - Sparse Identification of Nonlinear Dynamics
  - `1_SINDY_Hello_World.jl` - Basic SINDY example
  - `1_SINDY_Lorenz.jl` - SINDY with Lorenz system
  - `2_Universal_Differential_Equations_Dec2022.jl` - UDE examples
- **PINN/** - Physics-Informed Neural Networks
  - `0_BurgersEquation_MofLines.jl` - PINN for Burgers equation
- **DiffEqFlux/** - Differential equations with Flux
  - `Hello_NODE_GPU_NOTWORKING.jl` - GPU-based NODE (experimental)
- **UDE/** - Universal Differential Equations
- `NotWorking2022/` - Experimental SciML code

### 📊 **PDE_Ex/**
Partial Differential Equations and numerical methods.
- `0_HelloWorld_ModelingToolkit.jl` - Introduction to ModelingToolkit
- `1_MethodOfLines_Heat_Equation.jl` - Heat equation using method of lines
- `1_MethodOfLines_BurgersEquation.jl` - Burgers equation solution
- `DOCS_INFO_USEFUL/` - Documentation and reference materials
- `NotWorking/` - Experimental PDE implementations including:
  - Finite difference methods
  - Shallow water models
  - Wave equations
  - Poisson equations

### 📈 **Plots_Ex/**
Data visualization and plotting examples.
- `PlotsBasics.jl` - Comprehensive plotting examples with various backends
- `NotWorking/` - Experimental plotting code (Plotly, GR, Makie)

### 🌐 **Dash/**
Web application development with Dash.jl.
- `HelloDash.jl` - Basic Dash web application example
- `JuliaDash_OZ.txt` - Dash.jl usage notes

### ⚡ **CUDA/**
GPU computing with CUDA.jl.
- `Hello_CUDA_Notworking.jl` - CUDA basics (experimental)

### 🔧 **Optimization_Ex/**
Mathematical optimization examples.
- `OptimExamples.jl` - Various optimization problems and solvers

### 🌍 **GEO_Ex_Not_Tested_2022/**
Geospatial and geographic data processing.
- `GeoTiff_Examples.jl` - GeoTIFF file handling
- `MathBasics.jl` - Mathematical operations for geospatial data
- `ShapesAndPlots.jl` - Geometric shapes and plotting
- `GMT_Ex/` - Generic Mapping Tools examples
  - `GeneratePlots.jl` - GMT plotting utilities
  - `PlotMapsGMT.jl` - Map generation with GMT

### 📚 **NotOrganized_NotTested2022/**
Miscellaneous examples and experiments.
- `PlotsBasics.jl` - Basic plotting examples
- `NetCDF.jl` - NetCDF file handling
- `IO_basics.jl` - Input/output operations
- `AdvectionModel.jl` - Advection modeling
- `NavierStrokesPart1.jl` - Navier-Stokes equations
- `ODE.jl` - Basic ODE examples
- `ProfileExample.jl` - Code profiling
- **DataAssimilation/** - Data assimilation techniques
  - Kalman filters
  - Optimal interpolation
  - Gaussian processes

### 🛠️ **OZ_Tools_Not_Tested_2022/**
Custom tools and utilities.
- `DataAugmentation.jl` - Data augmentation techniques
- `ExampleData.jl` - Sample datasets
- **EOAS_OZLIB/** - Earth, Ocean, and Atmospheric Sciences library
  - `4DVisualizer.jl` - 4D data visualization
  - `GeoUtils.jl` - Geographic utilities
  - `NetCDF_Reader.jl` - NetCDF file reader
  - `RecurrentComputations.jl` - Recurrent computation patterns

### 🏛️ **Turin_Ex/**
Examples from Turin workshop/event.
- `BasicsNotWorking.jl` - Basic examples (experimental)

## 📋 Requirements

- Julia 1.6+ (recommended: latest stable version)
- Required packages are specified in each folder's `Project.toml`

## 🎯 Use Cases

This collection serves as a reference for:
- Learning Julia programming fundamentals
- Scientific computing and numerical methods
- Machine learning and neural networks
- Differential equations and mathematical modeling
- Data visualization and plotting
- Geospatial data processing
- Web application development
- GPU computing
- Optimization problems

## 📝 Notes

- Folders marked with "NotWorking" or "Not_Tested" contain experimental code
- Each example is designed to be self-contained and runnable independently
- Some examples may require specific Julia versions or package versions
- Check individual `Project.toml` files for exact dependencies

## 🤝 Contributing

Feel free to improve examples, add new ones, or fix issues in the experimental code sections.