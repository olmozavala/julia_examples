module OceananigansExamples

using CUDA
using Dash
using Dates
using LinearAlgebra
using NCDatasets
using Oceananigans
using Printf
using Statistics

include("catalog.jl")
include("runner.jl")
include("dashboard.jl")

export EXAMPLE_SPECS, RUNS_DIR, build_app, run_app, start_run!

end
