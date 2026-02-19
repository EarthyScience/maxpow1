 using Pkg # loads julia packagae manager
 Pkg.activate(".") # to keep the packages in the local environment; # ShareAdd is the shared environment
 Pkg.add("ModelingToolkit")

 
using ModelingToolkit # loads the package
isdefined(ModelingToolkit, Symbol("@mtkmodel")) # the safest check if macro exists in a package or not
using OrdinaryDiffEq
using Plots
# using RecursiveArrayTools # for VectorOfArray
# using DiffEqParamEstim
# -------------------------Model definition----------------------------------
# using ModelingToolkit: nothingasD, nothingast
using ModelingToolkit: t_nounits as t, D_nounits as D # imports specific symbols only from modeling toolkit
#t_nounits = symploic time variable with no units; D_nounits = symbolic derivative operator with no units
# @component FOL being ... end
@mtkmodel FOL begin # defines a modeling toolkit model name FOL
    @parameters begin
        m_A = 15 # IN THE PAPER, THIS IS BACK CALCULATED FROM STEADY STATE: gc/m2
        t_A =  12  # 12, 6 deacy rate / turnover rate yr^-1
        s_A = 0 # in half of month or 0.2 biomass will decay completely; yr^-1
        E = 0.4
        i = 400 # taking it at steady state
        l_e = 1
        k_Y = 12/12 #; maximum decomposition rate of the young pool yr^-1 ) (gone in 1 yr.)
        k_O = 1/100
        E_web = 0.4 # @ higher trophic levels e is more low bcz of movement, etc.
    end
    @variables begin # write them as function of independent variable
        Y(t) = 400 # pool: state variable dependent variables; initial condition at t=0
        O(t) = 6400 # pool: state variable THIS CHANGES VERY SLOWLY; Steady state is not achieved quixkly
        DOM(t) = 80 # pool: state variable
        A(t) = 4 # pool: state variable
        #RHS(t) #seperating right hand side calc, by introducing intermediate variable RHS
        d_Y(t) # flux: decomposition from young pool
        d_O(t) # flux: decomposition from old pool 
        u(t) # flux: uptake  @ steady state = sum of input decomposition fluxes 
        r_g(t) # flux: growth respiration
        r_m(t) # flux: maintenance respiration
        τ(t) # flux: microbial turnover backslash + ltish documentation (for greek letters code)
        l_A(t) # microbial limitation 
        D_A(t) # differential of A
        growth(t)
        r_web(t) # food web respiration
    end
    @equations begin
        D(Y) ~ (i - d_Y) # Pool 1: high-quality substrate; average input rate of 400gc/yr; not considering seasonal fluctuations of litter here
        D(O) ~ E_web*τ - d_O # Pool 2: low-quality substrate; added E_web term
        D(DOM) ~ d_Y + d_O - u # abstracted pool 3: assimilable @ steady stead this will be equal to 0
        D(A) ~ D_A
        D_A ~ u - r_g -r_m - τ # Pool 4: active microbial biomass
        d_Y ~ l_A * l_e * k_Y * Y # flux: decomposition of Y @ will this depend on Y present at any time t? or at t=0
        d_O ~ l_A * l_e * k_O * O # flux: decomposition of O
        u ~ d_Y + d_O # flux: uptake; Steady state expression; this is different from assim explicit model expression
        r_g ~ (1 - E) * u # flux: growth respiration; growth respiration is larger in our experiments
        r_m ~ s_A * A # flux: maintenance respiration ; flux of maintenance respiration is smaller than mcrobial turnover/ decay; building DNA, RNA
        τ ~ t_A * A # flux/pool? : microbial turnover
        l_A ~ A / (m_A + A) # MICROBIAL LIMITATION
        growth ~ u*E
        r_web ~ (1-E_web)*τ

    end
end

# 1) way turns our symbolic model into a numerical ODE problem and solves it  #main Julia package for numerically solving ODEs
@mtkcompile fol = FOL() # this creates the instance of my symbolic model fol is now a compiled ODE system, ready for numerical use

#--------------------------building base ODE problem------------------------------------------
prob = ODEProblem(fol, Dict(), (0.0, 10_000.0))# ODEProblem(system - compiled ODE system, u0 = []: initial conditions, tspan = time interval of simulation; p - parameters: here tau = 3.0)prob_2 = ODEProblem(fol, Dict(fol.k_Y => 4, fol.t_A => 16), (0.0,100.0))
#prob2 = ODEProblem(fol, Dict(flo.k_y => 4), (0.0, 100.0))
sol = solve(prob, reltol = 1e-9, abstol = 1e-9) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
# solution = sol(10_000, idxs = [fol.A, fol.D_A, fol.u*fol.E])

# plot(solution)#, u*E, τ, r_m)
# plot(sol.t => 10,000, sol( idxs=[fol.A, fol.D_A]) )
plot(sol, idxs=[:u, :D_A, :τ, :d_Y, :d_O, :r_m, :r_g, :l_A], xlabel="Time", ylabel="flux", title="State Variables over Time", legend =:outertopright)
# plot(sol.t,sol[fol.d_Y])
# plot(sol.t,sol[fol.A],xlabel="t", ylabel="A biomass", label = "A")
# plot(sol.t,sol[fol.D_A],xlabel="t", ylabel="dA/dt")
# plot(sol.t,sol[fol.u*fol.E],xlabel="t", ylabel="u uptake")
# plot(sol.t,sol[fol.τ],xlabel="t", ylabel="τ turnover")
# plot(sol.t,sol[fol.r_m],xlabel="t", ylabel="r_m maintenance")
#plot(sol.t,sol[fol.tau])
#plot(sol[fol.A])

#------------------------------------function to plot growth flux A vs t_A
function compute_steady_growth_t_A(t_A_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( fol.t_A => t_A_local))
    sol = solve(prob_local, maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
    sol[fol.τ][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #return (3-t_A_local)^2
end



function compute_steady_growth_E(E_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( fol.E => E_local))
    sol = solve(prob_local, 
    Rodas5(),reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
    sol[fol.growth][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #return (3-t_A_local)^2
end
function compute_steady_growth_E_web(E_web_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( fol.E_web => E_web_local))
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
    sol[fol.growth][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[fol.τ*E_web_local][end]
    #return (3-t_A_local)^2
end
#compute_steady_growth(3.0)

# Defining a vector variable
t_A = 0.01:0.05:20
k_Y = 1.0:0.4:20.0
E = 0.01:0.05:1.0
E_web = 0.01:0.05:0.6
pred_growth = compute_steady_growth_t_A.(t_A) 
pred_growth = compute_steady_growth_E.(E) # variable
pred_growth = compute_steady_growth_E_web.(E_web) # variable
pred_growth = compute_steady_growth_E(0.05) # variable: calling the function one time
#pred_growth = [compute_steady_growth.(t_A,k_Y) for  k_Y in k_Y_vals,t_A in t_A_vals] # variable
#plot( t_A,pred_growth, xlabel = "t", ylabel = "") # plot points 
plot( t_A,pred_growth, xlabel = "t_A", ylabel = "turnover", label = "turnover") # plot points 
plot( E,pred_growth, xlabel = "Efficiency", ylabel = "growth: uptake*E", label = "growth vs efficiency") # plot points 
plot( E_web,pred_growth)
#surface(t_A_vals, k_Y_vals, pred_growth,
        # xlabel="t_A", ylabel="k_Y", zlabel="Steady-state flux",
        # title="Steady-state flux as function of t_A and k_Y",
        # c=:viridis)# plot points 






#--------------------------generating symthetic data-------------------------------------------
# sol = solve(prob, Tsit5())
# ts = collect(range(0, stop = 100, length = 10))
# randomized = VectorOfArray([(sol(ts[i]) + 0.01randn(4)) for i in 1:length(t)])
# data = convert(Array, randomized)

#--------------------------building the cost function-----------------------------------
# cost_function = build_loss_objective(prob, Tsit5(), L2Loss(ts, data),
#     Optimization.AutoForwardDiff(),
#     maxiters = 10000, verbose = false)






# #--------------------------Parameter sweep for t_A--------------------------------------------
# tA_vals = 1.0:1.0:30.0      # Values of t_A to test
# Afinal = Float64[]           # Array to store results

# for val in tA_vals
#     newprob = remake(prob, p=[fol.t_A => val])   # change parameter
#     sol = solve(newprob, Tsit5(), saveat=100.0) # solve ODE: Tsit5 - Use this efficient 5th-order Runge-Kutta method; saveat - only stores the result at t = 100.0
#     push!(Afinal, sol[fol.A][end])              # store final A
# end

# # --------------------------Plot the results ----------------------------------------------
# plot(tA_vals, Afinal, xlabel="t_A", ylabel="Final A", lw=3,  title="Parameter Sweep for t_A")

# #---------------------------------------------------------------------------------------------------------------------
# sol = solve(prob) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
# #sol = solve(prob_2)

# plot(sol.t,sol[fol.d_Y])
#plot(sol[fol.A])


############## 
# see issue for clues:
# https://discourse.julialang.org/t/using-mtk-when-i-import-modelingtoolkit/133681/40