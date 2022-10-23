# Cool examples!!!!! https://tutorials.sciml.ai/html/models/01-classical_physics.html
############################################# Basics Calculus ####################################################
## Divergence (Vector to scalar)
∇⋅F = ∂Uₓ + ∂Vᵧ + ∂W𝕫   #→  For F(x,y,z) = (U(x,y,z), V(x,y,z), W(x,y,z))

## Laplacian (Scalar to scalar)
Δ²f = ∇⋅∇ = ∑ ∂²f/x²ᵢ  #→ For f(x₁, ...xₙ)

## Curl (Vector to vector) It uses n̂ as a parameter
curl F⋅n̂ = ∇ × F⋅n̂ = ((∂F𝕫/∂y - ∂F𝕪/∂𝕫), (∂F𝕫/∂y - ∂F𝕪/∂𝕫), (∂F𝕫/∂y - ∂F𝕪/∂𝕫))
# →  For F(x,y,z) = (U(x,y,z), V(x,y,z), W(x,y,z))

############################################# 2nd Order PDE classification ####################################################
Auₓₓ + Buₓ𝒚 + Cu𝒚𝒚 + I(x, y, u uₓ, u𝒚) = 0  # Where A, B and C are functions of x and y F(x,y)
# The classification depends on the values of f(x,y) = B² - 4AC
# f(x,y) > 0 → Hyperbolic
# f(x,y) = 0 → Parabolic
# f(x,y) < 0 → Elliptic


############################################# GENERAL EQUATIONS ####################################################
# ========= Poisson =========
Δu = f  #-->  ∇²u = f
Δu(x,y) = f(x,y)  #  In 2D

############################################# SPECIFIC EQUATIONS ####################################################
## ----- Advection equation
∂f/∂t + ∇ ⋅ (fu) = 0  → ∂f/∂t + (uᵢfₓ + uⱼf𝒚 + uₖf𝒛) = 0  
∂f/∂t = - ∇⋅(fu) → fₜ = (uᵢfₓ + uⱼf𝒚 + uₖf𝒛)

# -------- Wave equation
uₜₜ = c²uₓₓ

## ----- Heat/Diffusion (Poisson) -----
# Describes the flow of heat in a homogeneous and isotropic medium. 
# Homogeneous → Same stuff everywhere
# Isotropic → Same properties in all directions
uₜ = Δu
uₜ(t,x,y) = uₓₓ + u𝐲𝐲 # Heat
uₜ(t,x,y) = α(uₓₓ + u𝐲𝐲)   # Diffusion with α the diffusivity coefficient  

## ----- Convection/Diffusion or Advection/Diffusion equation
uₜ = ∇⋅(D∇u) - ∇⋅(𝒗u)+R  # D Diffusion coefficient, 𝒗 is the vector field, R describes sources or sinks
uₜ(t,x,y) = D(uₓₓ + u𝐲𝐲) - (𝒗ᵢuₓ + vⱼu𝒚) + R # With D constant