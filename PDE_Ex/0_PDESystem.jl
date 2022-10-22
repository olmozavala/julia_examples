using ModelingToolkit, MethodOfLines

@parameters x
@variables t u(..)
Dxx = Differential(x)^2
Dtt = Differential(t)^2
Dt = Differential(t)

## ------------- Wave equation -----------
C=1
eq  = Dtt(u(t,x)) ~ C^2*Dxx(u(t,x))

# Initial and boundary conditions
bcs = [u(t,0) ~ 0.,# for all t > 0
       u(t,1) ~ 0.,# for all t > 0
       u(0,x) ~ x*(1. - x), #for all 0 < x < 1
       Dt(u(0,x)) ~ 0. ] #for all  0 < x < 1]

# Space and time domains
domains = [t ∈ (0.0,1.0),
           x ∈ (0.0,1.0)]


@named pdesys = PDESystem(eq,bcs,domains,[t,x],[u(t,x)])

## Method of lines discretization
dx = 0.1
order = 2
discretization = MOLFiniteDifference([x=>dx],t, approx_order = 2)

## Convert the PDE problem into an ODE problem
prob = discretize(pdesys,discretization)

## Solve ODE problem
sol = solve(prob)

## Plot results
anim = @animate for i ∈ 1:length(sol.t)
plot(sol.u[i], label = "wave", ylims =[-1, 1])
end every 5

gif(anim, "1Dwave.gif", fps = 10)