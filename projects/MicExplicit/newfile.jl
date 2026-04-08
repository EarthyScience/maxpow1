 using Pkg # loads julia packagae manager
 Pkg.activate(".") # to keep the packages in the local environment; # ShareAdd is the shared environment
 Pkg.add("ModelingToolkit")

 
using ModelingToolkit # loads the package
isdefined(ModelingToolkit, Symbol("@mtkmodel")) # the safest check if macro exists in a package or not
using OrdinaryDiffEq
using Plots
using SteadyStateDiffEq
using DifferentialEquations
# using RecursiveArrayTools # for VectorOfArray
# using DiffEqParamEstim

# -------------------------Model definition----------------------------------
# using ModelingToolkit: nothingasD, nothingast
using ModelingToolkit: t_nounits as t, D_nounits as D # imports specific symbols only from modeling toolkit
#t_nounits = symploic time variable with no units; D_nounits = symbolic derivative operator with no units
# @component ofxd being ... end
@mtkmodel OFXD begin # defines a modeling toolkit model name ofxd
    @parameters begin
        m_A = 15 # IN THE PAPER, THIS IS BACK CALCULATED FROM STEADY STATE: gc/m2
        t_A =  12  # 12, 6 deacy rate / turnover rate yr^-1
        t_E = 5 
        s_A = 0 # in half of month or 0.2 biomass will decay completely; yr^-1
        ϵ = 0.5
        #E_c = 0.4
        #E_e = 0.5
        i = 400 # taking it at steady state
        l_e = 1
        k_Y = 12/12 #; maximum decomposition rate of the young pool yr^-1 ) (gone in 1 yr.)
        k_O = 1/100
        E_web = 0.4 # @ higher trophic levels e is more low bcz of movement, etc.
        a = 0.2
        k_abio_O = k_O/1000 # abiotic loss constants
        k_abio_Y = k_Y/1000 # abiotic loss constants
    end
    @variables begin # write them as function of independent variable
        Y(t) = 400 # pool: state variable dependent variables; initial condition at t=0
        O(t) = 6400 # pool: state variable THIS CHANGES VERY SLOWLY; Steady state is not achieved quixkly
        DOM(t) = 80 # pool: state variable
        A(t) = 4 # pool: state variable
        #E(t) = 0.3 # pool: state variable enzymes pool
        #RHS(t) #seperating right hand side calc, by introducing intermediate variable RHS
        D_A(t) # differential of A
        #D_E(t) # differential of E
        d_Y(t) # flux: decomposition from young pool
        d_O(t) # flux: decomposition from old pool 
        dt_O(t) # O pool changing with time
        u(t) # flux: uptake  @ steady state = sum of input decomposition fluxes 
        #u_E(t)
        r_g(t) # flux: growth respiration
        r_m(t) # flux: maintenance respiration
        #r_E(t) # flux: enzyme respiration
        #τ(t) # flux: microbial turnover backslash + ltish documentation (for greek letters code)
        τ_A(t)
        τ_E(t)
        l_A(t) # microbial limitation 
        growth_A(t)
        A1(t)
        growth(t)
        growth_E(t)
        r_web(t) # food web respiration
        r_total(t) # total respiration: check   
        abio_Y(t) # abiotic flux from y pool
        abio_O(t) # abiotic flux from o pool
        abio_sum(t)
    end
    @equations begin
        D(Y) ~ (i - d_Y - abio_Y) # Pool 1: high-quality substrate; average input rate of 400gc/yr; not considering seasonal fluctuations of litter here
        dt_O ~ E_web*τ_A + τ_E - d_O - abio_O# Pool 2: low-quality substrate; added E_web term
        D(O) ~ 0 # it is not in seady state but we ignore the change
        D(DOM) ~ d_Y + d_O - u # abstracted pool 3: assimilable @ steady stead this will be equal to 0
        D(A) ~ D_A
        #D(E) ~ D_E     # enzyme pool 
        #D_A ~ u  - r_g -r_m - τ # Pool 4: active microbial biomass
        D_A ~ u - r_g -r_m - growth_E - τ_A # u-(1-e_c)*u
        #D_E ~ u_E - τ_E 
        d_Y ~ l_A * l_e * k_Y * Y # flux: decomposition of Y @ will this depend on Y present at any time t? or at t=0
        d_O ~ l_A * l_e * k_O * O # flux: decomposition of O
        u ~ d_Y + d_O # flux: uptake; Steady state expression; this is different from assim explicit model expression
        #u_E = (a*A)/ϵ
        r_g ~ (1 - ϵ) * (u) # flux: growth respiration; growth respiration is larger in our experiments
        #r_g ~ (1-ϵ)*(u - u_E)
        # #r_E ~ (1-ϵ)*(u_E)
        r_m ~ s_A * A # flux: maintenance respiration ; flux of maintenance respiration is smaller than mcrobial turnover/ decay; building DNA, RNA
        #τ ~ t_A * A # flux/pool? : microbial turnover
        τ_A ~ t_A * A
        τ_E ~  growth_E# assumption E is at steady state 
        l_A ~ (a*A) / (m_A + a*A) # MICROBIAL LIMITATION
        A1 ~ a*A
        growth ~ u*ϵ
        growth_A ~ u * ϵ - growth_E
        growth_E ~ a*A # Assumption: Growth (production of new enezymes) is proportional to A (microbial biomass). 
        r_web ~ (1 - E_web)*τ_A
        r_total ~ r_g + r_m + r_web + abio_sum
        abio_sum ~ abio_O +abio_Y
        abio_O ~  k_abio_O*O
        abio_Y ~ k_abio_Y*Y
    end
end

# 1) way turns our symbolic model into a numerical ODE problem and solves it  #main Julia package for numerically solving ODEs
@mtkcompile ofxd = OFXD() # this creates the instance of my symbolic model ofxd is now a compiled ODE system, ready for numerical use

#--------------------------building base ODE problem------------------------------------------
#Dict(ofxd.a => 0.141)
# ODEProblem(system - compiled ODE system, u0 = []: initial conditions, tspan = time interval of simulation; p - parameters: here tau = 3.0)prob_2 = ODEProblem(ofxd, Dict(ofxd.k_Y => 4, ofxd.t_A => 16), (0.0,100.0))
prob = ODEProblem(ofxd, Dict(), (0.0, 10_000.0))
sol = solve(prob, reltol = 1e-9, abstol = 1e-9) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
plot(sol)
u0 = sol.u[end]
sol[ofxd.dt_O][end]


plot(sol, idxs =[:growth_A])
plot(sol, idxs=[:u, :D_A, :τ_A, :d_Y, :d_O, :r_m, :r_g, :l_A, :r_total], xlabel="Time", ylabel="flux", title="State Variables over Time", legend =:outertopright)
sol[ofxd.r_total][end]
sol[ofxd.abio_sum][end]
# plot(sol.t,sol[ofxd.d_Y])
# plot(sol.t,sol[ofxd.A],xlabel="t", ylabel="A biomass", label = "A")
# plot(sol.t,sol[ofxd.D_A],xlabel="t", ylabel="dA/dt")
# plot(sol.t,sol[ofxd.u*ofxd.E],xlabel="t", ylabel="u uptake")
# plot(sol.t,sol[ofxd.τ],xlabel="t", ylabel="τ turnover")
# plot(sol.t,sol[ofxd.r_m],xlabel="t", ylabel="r_m maintenance")
#plot(sol.t,sol[ofxd.tau])
plot(sol[ofxd.A])

#------------------------------------function to plot growth flux A vs t_A
function compute_steady_growth_t_A(t_A_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( ofxd.t_A => t_A_local))
    sol = solve(prob_local, maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
    #sol = solve(prob_2)
    sol[ofxd.τ_A][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #return (3-t_A_local)^2
end



function compute_steady_growth_E(E_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( ofxd.E => E_local))
    sol = solve(prob_local, Rodas5(),reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
    #sol = solve(prob_2)
    sol[ofxd.growth][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #return (3-t_A_local)^2
end
function compute_steady_growth_E_web(E_web_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local = remake(prob, p = Dict( ofxd.E_web => E_web_local))
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
    #sol = solve(prob_2)
    sol[ofxd.growth][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[ofxd.τ*E_web_local][end]
    #return (3-t_A_local)^2
end
#compute_steady_growth(3.0)

function compute_steady_growth_A(E_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local0 = remake(prob, p = Dict( ofxd.ϵ => E_local))
    prob_local = remake(prob_local0; u0)
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
 #sol = solve(prob_2)
    sol[ofxd.growth_A][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[ofxd.τ*E_web_local][end]
    #return (3-t_A_local)^2
end

function compute_steady_growth_A(a_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local0 = remake(prob, p = Dict( ofxd.a => a_local))
    prob_local = remake(prob_local0; u0)
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
 #sol = solve(prob_2)
    sol[ofxd.growth_A][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[ofxd.τ*E_web_local][end]
    #return (3-t_A_local)^2
end

function compute_steady_growth(a_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local0 = remake(prob, p = Dict( ofxd.a => a_local))
    prob_local = remake(prob_local0; u0)
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
    sol[ofxd.growth][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[ofxd.τ*E_web_local][end]
    #return (3-t_A_local)^2
end

function compute_steady_growth_1(a_local::Float64)
    # remake the ODE problem with multiple parameters
    prob_local0 = remake(prob, p = Dict( ofxd.a => a_local))
    prob_local = remake(prob_local0; u0)
    sol = solve(prob_local, 
   # Rodas5(),
    reltol = 1e-9, abstol = 1e-9,
      maxiters = Int(1e7), saveat=[0,10_000]) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
#sol = solve(prob_2)
    sol[ofxd.A1][end] # [end]: this takes the last valuethe uptake flux will be equal to the turnover flux. what we are maximising is the uptake flux (bcz in steady state = turnover so we are computing turnover)
    #sol[ofxd.τ*E_web_local][end]
    #return (3-t_A_local)^2
end


# Defining the analytical function parameters
 #E_web = 0.4
#  t_A = 12.0
#  a = 0.2
#  m_A = 15.0
#  l_e = 1.0
#  k_abio_O = 1/10^5

# Defining a vector parameters
t_A = 0.01:0.05:20
k_Y = 1.0:0.4:20.0
E = 0.01:0.05:1.0
ϵ = 0.01:0.05:1.0
E_web = 0.01:0.05:0.6
#a = 0.001:0.02:0.07 # can assosciate range with no. of values
a = range(0.001, 0.07, length = 50)
a[7]
#a = 0.05
#
# function analytical(; O::Float64,E_web::Float64, k_O::Float64, t_A::Float64, a::Float64, m_A::Float64, l_e::Float64, k_abio_O::Float64) # semicolon means named arguments
#function analytical(; O::Float64, i::Float64, ϵ::Float64, k_Y::Float64, k_O::Float64, t_A::Float64, a::Float64, m_A::Float64, s_A::Float64, k_abio_Y::Float64) # semicolon means named arguments
pdef = (;
    O = convert(Float64, defaults(ofxd)[ofxd.O]),
    i = convert(Float64, defaults(ofxd)[ofxd.i]),
    ϵ = convert(Float64, defaults(ofxd)[ofxd.ϵ]),
    k_Y = convert(Float64, defaults(ofxd)[ofxd.k_Y]),
    k_O = convert(Float64, defaults(ofxd)[ofxd.k_O]),
    t_A = convert(Float64, defaults(ofxd)[ofxd.t_A]),
    a = convert(Float64, defaults(ofxd)[ofxd.a]),
    m_A = convert(Float64, defaults(ofxd)[ofxd.m_A]),
    s_A = convert(Float64, defaults(ofxd)[ofxd.s_A]),
    k_abio_Y = 1/10^3 #convert(Float64, defaults(ofxd)[ofxd.k_abio_O]),
)
#convert.(Float64, pdef)   
#include("function.jl")
include("Function_2.jl")
#pred_steady_1 = analytical(;pdef...)
pred_steady_2 = analytical(;pdef...)
if(pred_steady_2 > 0) 
    Young = i/((a*pred_steady_2)/(m_A+a*pred_steady_2) + k_abio_Y)
end
#pred_steady_1 = analytical(E_web, t_A, a, m_A, l_e, k_abio_O)
# pred_steady_2 = analytical(i, ϵ,k_Y,k_O, t_A, a, m_A, s_A, k_abio_Y)
# pred_growth = compute_steady_growth_t_A.(t_A) 
pred_growth = compute_steady_growth_A.(a)
pred_growth = compute_steady_growth.(a)
# pred_growth_1 = compute_steady_growth_1.(a)
# pred_growth = compute_steady_growth_t_A.(t_A) 
# pred_growth = compute_steady_growth_A.(a1)
# pred_growth = compute_steady_growth_A.(a2) 
# pred_growth = compute_steady_growth_A.(ϵ) # variable
# pred_growth = compute_steady_growth_E_web.(E_web) # variable
# pred_growth = compute_steady_growth_E(0.05) # variable: calling the function one time
#pred_growth = [compute_steady_growth.(t_A,k_Y) for  k_Y in k_Y_vals,t_A in t_A_vals] # variable
#plot( t_A,pred_growth, xlabel = "t", ylabel = "") # plot points 
# plot( t_A,pred_growth, xlabel = "t_A", ylabel = "turnover", label = "turnover") # plot points 
# plot( E,pred_growth, xlabel = "Efficiency", ylabel = "growth: uptake*E", label = "growth vs efficiency") # plot points 
# plot( E_web,pred_growth)
# plot( ϵ,pred_growth)
# plot( a,pred_growth)
#  plot(a, pred_growth, seriestype =:scatter)
#plot(a[4:end], pred_growth[4:end], seriestype =:scatter, xlabel = "yield coefficient", ylabel = "Microbial growth", label = "growth vs enzyme yield coefficient")
plot(a[4:end], pred_growth[4:end], linewidth=3, guidefontsize=14, legendfontsize=14,tickfontsize=14, xlabel = "specific enzyme allocation [yr^-1]", ylabel = "Microbial growth [gcm^-2yr^-1]" ,
            grid =false, dpi = 300, title = "Microbial growth vs Specific enzyme allocation", legend = false, titlefont = font(14))
# imax = argmax(pred_growth)     # index of maximum y
# x_max = a[imax]      # corresponding x value
# y_max = pred_growth[imax]      # maximum y value
p = plot(a[4:end], pred_growth[4:end], linewidth=3, guidefontsize=14, legendfontsize=14,tickfontsize=14, xlabel = "specific enzyme allocation [yr^-1]", ylabel = "Microbial growth [gcm^-2yr^-1]" ,
            grid =false, dpi = 300, title = "Microbial growth vs Specific enzyme allocation", legend = false, titlefont = font(14))
#plot!(a[4:end], pred_growth_1[4:end])
plot!(twinx(), a[4:end], pred_growth_1[4:end], label="Enzyme investment", ylabel="Enzyme investment", color=:red)
dA/dt = ϵ*u(a) - a*A - τ_A
#surface(t_A_vals, k_Y_vals, pred_growth,
        # xlabel="t_A", ylabel="k_Y", zlabel="Steady-state flux",
        # title="Steady-state flux as function of t_A and k_Y",
        # c=:viridis)# plot points 
# plot(x, y,
#     linewidth=3,
#     label="Entropy Production",
#     xlabel="Time (s)",
#     ylabel="Ṡ_gen",
#     title="Entropy Production vs Time",
#     legend=:topright,
#     framestyle=:box,
#     grid=false,
#     dpi=300
# )





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
#     newprob = remake(prob, p=[ofxd.t_A => val])   # change parameter
#     sol = solve(newprob, Tsit5(), saveat=100.0) # solve ODE: Tsit5 - Use this efficient 5th-order Runge-Kutta method; saveat - only stores the result at t = 100.0
#     push!(Afinal, sol[ofxd.A][end])              # store final A
# end

# # --------------------------Plot the results ----------------------------------------------
# plot(tA_vals, Afinal, xlabel="t_A", ylabel="Final A", lw=3,  title="Parameter Sweep for t_A")

# #---------------------------------------------------------------------------------------------------------------------
# sol = solve(prob) # sol is an ODESolution object; chooses a default ODE solver sol conatins sol.t; sol.u; sol(t) - interpolated solution
# #sol = solve(prob_2)

# plot(sol.t,sol[ofxd.d_Y])
#plot(sol[ofxd.A])


############## 
# see issue for clues:
# https://discourse.julialang.org/t/using-mtk-when-i-import-modelingtoolkit/133681/40