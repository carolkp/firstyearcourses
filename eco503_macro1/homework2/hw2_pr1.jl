# Homework 2

#Problem 1

# Functions

# Calculate Output
function output(k, A, alpha)
	return A*k^(alpha)
end

# Calculate derivative of output function
function output_derivative(k, A, alpha)
	return alpha*A*k^(alpha - 1)
end

# Find capital value in steady state
function find_steady_state(A, alpha, rho, delta)
	aux_power = 1/(1 - alpha)
	aux_fraction = alpha*A / (rho + delta)
	return aux_fraction^aux_power
end

# Calculate next value of capital with discrete approximation of derivative
function next_capital(previous_k, previous_c, step_size, A, alpha, delta)
	aux_product = output(previous_k, A, alpha) - delta*previous_k - previous_c
	return previous_k + step_size*aux_product
end

# Calculate next value of consumption with discrete approximation of derivative
function next_consumption(previous_k, previous_c, step_size, A, alpha, delta, sigma, rho)
	aux_product = rho + delta - output_derivative(previous_k, A, alpha)
	return previous_c*(1 - (step_size/sigma)*aux_product)
end

# Find best initial consumption
function calculate_path(number_of_steps, initial_k, initial_c)
	# Initialize sequence of capital and consumption for current iteration
	k_evolution = zeros(number_of_steps + 1, 1)
	k_evolution[1] = initial_k

	c_evolution = zeros(number_of_steps + 1, 1)
	c_evolution[1] = initial_c
	# Save sequence of errors
	error_evolution = ones(number_of_steps, 1)

	for index in 1:(number_of_steps)
		previous_k = k_evolution[index]
		previous_c = c_evolution[index]

		k_evolution[index + 1] = next_capital(previous_k, previous_c, step_size, A, alpha, delta)
		c_evolution[index + 1] = next_consumption(previous_k, previous_c, step_size, A, alpha, delta, sigma, rho)
		
		#Non-negativity constraints
		k_evolution[index + 1] = max(k_evolution[index+1], 0.1*tolerance)
		c_evolution[index + 1] = max(c_evolution[index+1], 0.1*tolerance)
		# Feasibility of consumption
		c_evolution[index + 1] = min(c_evolution[index + 1], output(k_evolution[index + 1], A, alpha) - delta * k_evolution[index + 1] + (1/step_size)*k_evolution[index + 1])
		error_evolution[index] = abs(k_steady_state - k_evolution[index + 1])
	end
	return k_evolution, c_evolution, error_evolution
end

function find_best_initial_consumption(tolerance, max_iterations, initial_k, initial_c)
	# Initialize iteration counter and variable to store error
	iteration_number = 1
	_error = 1
	# We'll test several different values for initial consumption
	# Leave loop if after 500 iterations we're close to steady state 
	#(or if we reach the maximum number of iterations without converging)
	
	while ((iteration_number <= max_iterations) && (_error > tolerance))
	    println("\nIteration number: ", iteration_number)

		k_evolution, c_evolution, error_evolution = calculate_path(number_of_steps, initial_k, initial_c)
		_error = error_evolution[end]

		# Update initial guess for c according to last capital value (compared to steady state value)
		if (k_evolution[end] > k_steady_state)
			initial_c = (1 + 0.5/iteration_number)*initial_c
		else
			initial_c = (1  - 0.5/iteration_number)*initial_c
		end

		println("Error Value: ", _error)
		iteration_number += 1

		if ((_error <= tolerance) || (iteration_number > max_iterations))
			global k_path = k_evolution
			global c_path = c_evolution
			global error_path = error_evolution
		end
	end

	return k_path, c_path, error_path
end

# Given Parameters:
A = 1
alpha = 1/3
beta = 0.96
rho = 1/beta - 1
delta = 0.08
sigma = 1.01
step_size = 0.05
number_of_steps = 600
tolerance = 0.001
max_iterations = 1000

# Calculate initial capital level from initial steady state
k_steady_state = find_steady_state(A, alpha, rho, delta)
initial_k = 0.5*k_steady_state
#initial_c = output(initial_k, A, alpha) - delta*initial_k # just an initial guess
initial_c = 1

# Find initial consumption for saddle path with predefined initial k
k_path_best, c_path_best, error_path_best = find_best_initial_consumption(tolerance, max_iterations, initial_k, initial_c)

# Print Results
println()
println("RESULTS:")
println("\nBest initial consumption value: ", c_path_best[1])
println("Approximate Consumption in steady state: ", c_path_best[end])
println("Approximate Capital in steady state: ", k_path_best[end])
println("Steady State Capital: ", k_steady_state)
println("Final error value: ", error_path_best[end])



# Part b
c_best = c_path_best[1]
c_above = 1.5*c_best
c_below = 0.5*c_best
time_grid = range(0, step = step_size, length = number_of_steps+1)


k_path_above, c_path_above, _ = calculate_path(number_of_steps, initial_k, c_above)
k_path_below, c_path_below, _ = calculate_path(number_of_steps, initial_k, c_below)


all_c_paths = hcat(c_path_best, c_path_below, c_path_above)
all_k_paths = hcat(k_path_best, k_path_below, k_path_above)

using Plots
p_k = plot(time_grid, all_k_paths, lw = 1.5, title="Time Series for capital", legend = false)
p_c = plot(time_grid, all_c_paths, lw = 1.5, title="Time Series for consumption", label = ["Best Consumption" "Low Consumption" "High Consumption"], legend = :outerbottom, xlabel = "Time")
# display(plot(p_k, p_c, layout = grid(2,1, heights=[0.3, 0.7])))



p = plot(p_k, p_c, layout = grid(2,1, heights=[0.3, 0.7]), size = (1000, 800))
savefig(p, "1b_timeseries.png")

p2 = plot([k_path_best k_path_above k_path_below], [c_path_best c_path_above c_path_below], lw = 1.5, label = ["Best Consumption" "High Consumption" "Low Consumption"])
savefig(p2, "1b_diagram")