#Name: Marcelo Barbosa Ferreira
#Name: Carolina Kowalski Piazza
#Macro 1 Part 2 Homework 4


using Plots
using StatsBase
using Optim
using Statistics
using JLD

# Inputs: mu, rho, sig of an AR(1) process y(t) = rho y(t-1) + e(t)
# mu is the unconditional mean of the process
# rho is the AR coefficient
# sig is the standard deviation of e(t)

# Summary: The function discretizes the AR(1) process into an n equally
# spaced state Markov chain with transition probabilities given by the
# Rouwenhorst (1995) method. States are in the set [mu - nu , mu + nu]

# Outputs: Lambda, P
# Lambda is a nx1 vector of equally spaced states centered around mu
# P is the Markov transition matrix
function rouwenhorst(mu, rho, sig, n)

    nu = sqrt( ((n-1)/(1-rho^2)))*sig

    Lambda = collect(range(mu-nu, stop = mu+nu, length = n))

    p = (1+rho)/2
    q = p

    P1 = [ p 1-p; 1-q q ]

    if n == 2
        P = P1
    else
        for ii = 3:n
            zcol = zeros(ii-1,1)
            zrow = zcol'

            A = [ P1 zcol ; zrow 0 ]
            B = [ zcol P1 ; 0 zrow ]
            C = [ zrow 0 ; P1 zcol ]
            D = [ 0 zrow ; zcol P1 ]

            P1 = p*A + (1-p)*B + (1-q)*C + q*D
            P1[2:end-1,:] = P1[2:end-1,:]/2
        end
        P = P1
    end

    return Lambda, P
end



function utility(current_k, next_k, shock, labor)
	log(max(exp(shock)*(current_k^alpha)*labor^(1-alpha) + (1-delta)*current_k - next_k - theta*(labor^(1 + gamma))/(1+gamma),0))
end

function consumption(current_k, next_k, shock, labor)
	exp(shock)*(current_k^alpha)*labor^(1-alpha) + (1-delta)*current_k - next_k
end

function output(current_k, shock, labor)
	exp(shock)*(current_k^alpha)*labor^(1-alpha)
end

function find_labor(current_k, shock)
    (((1-alpha)/theta)*exp(shock)*current_k^alpha)^(1/(gamma + alpha))
end

function find_wage(current_k, labor, shock)
    exp(shock)*(1-alpha)*(current_k^(alpha))*(labor^(-alpha))
end


# Model Parameters
import Random
Random.seed!(1234)
alpha=0.4
beta=0.98
delta=0.01
xi = 0.96
sigma = 0.008
n_ss = 1/3
gamma = 1

# Numerical Parameters
tol = 0.00001
grid_size = 201
grid_tfp = 21



# Finding theta
ratio_kn = ( (1/alpha)*(1/beta - 1 + delta) )^(1/(alpha - 1))
ratio_cn = ratio_kn^alpha - delta*ratio_kn

theta = ((1-alpha)/n_ss^gamma)*ratio_kn^alpha

# Steady State Values
k_ss = ratio_kn * n_ss
c_ss = ratio_cn * n_ss



# Create grids for capital and possible shocks
gridz, transition_matrix = rouwenhorst(0,xi,sigma, grid_tfp)

kmin = 0.5*k_ss
kmax = 1.5*k_ss
gridk = LinRange(kmin,kmax,grid_size)


#Initial value of Value Function
V = zeros(grid_size,grid_tfp)
Vaux = zeros(grid_size,grid_tfp,grid_size)

#Policy rule vector
policyrule = zeros(grid_size,grid_tfp)



optimal_labor = zeros(grid_size,grid_tfp,grid_size)

println("\n Start calculating Value Function")

#Value function iteration
resid = 1
iter = 0

while resid>tol && iter<1000
    global iter=iter+1
    if iter%10 == 0
    	println(iter)
    end
    for i=1:grid_size
        for j=1:grid_tfp
            for l=1:grid_size
                optimal_labor[i,j,l] = find_labor(gridk[i], gridz[j])
                global Vaux[i,j,l]= utility(gridk[i], gridk[l], gridz[j], optimal_labor[i,j,l])  + beta*sum(transition_matrix[j,:].*V[l,:])
            end
        end
    end
    global resid=0
    for i=1:grid_size
        for j=1:grid_tfp
            value, ind = findmax(Vaux[i,j,:])
            global resid=max(resid,abs(value-V[i,j]))
            global V[i,j]=value
            global policyrule[i,j]=ind
        end
    end
end
println(iter)

policyrule = round.(Int,policyrule)


save("hw4_q3-value_function.jld", "value", V)
save("hw4_q3-policy_rule.jld", "policy", policyrule)

z_values_idx = round.(Int,[1,11,21])
labor_policy_function = zeros(grid_size, 3)
policy_function = zeros(grid_size, 3)
for i = 1:grid_size
    for j = 1:3
        labor_policy_function[i,j] = optimal_labor[i, z_values_idx[j], policyrule[i,z_values_idx[j]]]
        policy_function[i,j] = gridk[policyrule[i, z_values_idx[j]]]
    end
end

plot(gridk, labor_policy_function, title="Labor Policy Function", label=["z = -0.128" "z = 0" "z = 0.128"], lw=2, xlabel = "capital", ylabel="labor", legend= :bottomright)
png("hw4_q3-labor_policy")

plot(gridk, policy_function, title="Policy Function", label=["z = -0.128" "z = 0" "z = 0.128"], lw=2, xlabel = "current capital", ylabel="next capital", legend= :bottomright)
png("hw4_q3-capital_policy")

##############################################################################################

T = 500
shocks_idx = zeros(T)
z0 = findfirst(x->x==0, gridz)
shocks_idx[1] = z0
shocks_idx = round.(Int,shocks_idx)

for t = 2:T
    shocks_idx[t] = round(Int, sample(1:grid_tfp, Weights(transition_matrix[shocks_idx[t-1],:])) )
end

#Transition paths

k_index_path=ones(T)

kpath=k_ss*ones(T)
k0 = findfirst(x -> x >= k_ss, gridk)
k_index_path[1] = k0
k_index_path = round.(Int,k_index_path)
policyrule = round.(Int,policyrule)

cons_path = c_ss*ones(T)
output_path = (k_ss^(alpha))*(n_ss^(1-alpha))*ones(T)
labor_path = n_ss*ones(T)

shock_path = zeros(T)

for i=2:T
    shock_path[i] = gridz[shocks_idx[i]]
    k_index_path[i] = policyrule[k_index_path[i-1],shocks_idx[i-1]]
    kpath[i]=gridk[k_index_path[i]]
    labor_path[i] = optimal_labor[k_index_path[i],shocks_idx[i],policyrule[k_index_path[i],shocks_idx[i]]]
    cons_path[i] = consumption(kpath[i], gridk[policyrule[k_index_path[i],shocks_idx[i]]], gridz[shocks_idx[i]], labor_path[i])
    output_path[i] = output(kpath[i], gridz[shocks_idx[i]], labor_path[i])
end

plot(1:T, shock_path, title="Shock Path", lw=2, legend=false)
png("hw4_q3-shock_path")

plot(1:T, kpath,title="Capital Path",lw=2, legend=false)
png("hw4_q3-capital_path")

plot(1:T, output_path,title="Output path",lw=2, legend=false)
png("hw4_q3-output_path")

plot(1:T, labor_path,title="Labor path",lw=2, legend=false)
png("hw4_q3-labor_path")

plot(1:T, cons_path,title="Consumption path",lw=2, legend=false)
png("hw4_q3-consumption_path")
# plot(1:T, [shock_path, labor_path])

desc_statistics = zeros(4,3)
data = hcat(kpath, cons_path, labor_path, output_path)

for i = 1:4
    desc_statistics[i,1] = std(data[:,i])
    desc_statistics[i,2] = autocor(data[:,i], [1])[1]
    desc_statistics[i,3] = cor(data[:,i], data[:,4])
end

# Item 5

persistent_shock = zeros(T)

for t = 1:T
    persistent_shock[t] = findfirst(x -> x >= sigma*xi^(t-1), gridz)
end

persistent_shock = round.(Int,persistent_shock)

#Transition paths

k_index_path_pers=ones(T)
kpath_pers=k_ss*ones(T)

k0_pers = findfirst(x -> x >= k_ss, gridk)
k_index_path_pers[1] = k0_pers
k_index_path_pers = round.(Int,k_index_path_pers)

cons_path_pers = c_ss*ones(T)
output_path_pers = (k_ss^(alpha))*(n_ss^(1-alpha))*ones(T)
labor_path_pers = n_ss*ones(T)

shock_path_pers = zeros(T)
shock_path_pers[1] = gridz[persistent_shock[1]]

wage_pers = find_wage(k_ss, n_ss, shock_path_pers[1])*ones(T)

for i=2:T
    shock_path_pers[i] = gridz[persistent_shock[i]]
    k_index_path_pers[i] = policyrule[k_index_path_pers[i-1],persistent_shock[i-1]]
    kpath_pers[i]=gridk[k_index_path_pers[i]]
    labor_path_pers[i] = optimal_labor[k_index_path_pers[i],persistent_shock[i],policyrule[k_index_path_pers[i],persistent_shock[i]]]
    cons_path_pers[i] = consumption(kpath_pers[i], gridk[policyrule[k_index_path_pers[i],persistent_shock[i]]], gridz[persistent_shock[i]], labor_path_pers[i])
    output_path_pers[i] = output(kpath_pers[i], gridz[persistent_shock[i]], labor_path_pers[i])
    wage_pers[i] = find_wage(kpath_pers[i], labor_path_pers[i], shock_path_pers[i])
end

plot(1:T, shock_path_pers, title="Shock Path", lw=2, legend=false)
png("hw4_q3-shock_path_pers")

plot(1:T, kpath_pers,title="Capital Path",lw=2, legend=false)
png("hw4_q3-capital_path_pers")

plot(1:T, output_path_pers,title="Output path",lw=2, legend=false)
png("hw4_q3-output_path_pers")

plot(1:T, labor_path_pers,title="Labor path",lw=2, legend=false)
png("hw4_q3-labor_path_pers")

plot(1:T, cons_path_pers,title="Consumption path",lw=2, legend=false)
png("hw4_q3-consumption_path_pers")

plot(1:T, wage_pers,title="Wage path",lw=2, legend=false)
png("hw4_q3-wage_path_pers")

#