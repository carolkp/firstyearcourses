# Homework 2

# Question 3
# This question has you compare the optimal solution to the path produced by the Solow model. In
# particular, assume the same calibration as in (1), but let k_0 = k*/5. 
# Solve for the steady state corresponding to the optimal growth problem. 
# Solve for the fraction of output that is devoted to investment in the steady state, and denote this ratio by s. 
# Solve for both the optimal transition path and the transition path that is implied by the Solow model assuming a savings rate equal to the value
# of s that you just calculated. 
# Plot the two series for k_t on a single plot. Describe the differences. 
# Try to explain why the optimal solution deviates from the Solow solution in the manner that it does. You
# may find it useful to plot the two series for consumption.

# Calculate Output
function output(k)
	return A*(k^alpha)
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
	aux_product = output(previous_k) - delta*previous_k - previous_c
	return previous_k + step_size*aux_product
end

function next_capital_solow(previous_k, previous_c, step_size, A, alpha, delta, savings)
	aux_product = (savings*output(previous_k)) - (delta*previous_k)
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
		c_evolution[index + 1] = min(c_evolution[index + 1], output(k_evolution[index + 1]) - delta * k_evolution[index + 1] + (1/step_size)*k_evolution[index + 1])
		error_evolution[index] = abs(k_steady_state - k_evolution[index + 1])
	end
	return k_evolution, c_evolution, error_evolution
end

function calculate_path_solow(number_of_steps, initial_k, initial_c, savings)
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

		k_evolution[index + 1] = next_capital_solow(previous_k, previous_c, step_size, A, alpha, delta, savings)
		k_evolution[index + 1] = max(k_evolution[index+1], 0.1*tolerance)

		c_evolution[index + 1] = (1 - savings)*output(k_evolution[index + 1])
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
	    #println("\nIteration number: ", iteration_number)

		k_evolution, c_evolution, error_evolution = calculate_path(number_of_steps, initial_k, initial_c)
		_error = error_evolution[end]

		# Update initial guess for c according to last capital value (compared to steady state value)
		if (k_evolution[end] > k_steady_state)
			initial_c = (1 + 0.5/iteration_number)*initial_c
		else
			initial_c = (1  - 0.5/iteration_number)*initial_c
		end

		#println("Error Value: ", _error)
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
number_of_steps = 1000
tolerance = 0.001
max_iterations = 1000

# Calculate initial capital level from initial steady state
k_steady_state = find_steady_state(A, alpha, rho, delta);
initial_k = 0.2*k_steady_state;
initial_c_best = 1;

# Find optimal path from Neoclassical Growth Model (NGM)
k_path_best, c_path_best, error_path_best = find_best_initial_consumption(tolerance, max_iterations, initial_k, initial_c_best);

#Savings for Solow Model
savings = delta*k_steady_state/output(k_steady_state);
initial_c_solow = (1 - savings) * output(initial_k);

# Calculate path for Solow
k_path_solow, c_path_solow, _ = calculate_path_solow(number_of_steps, initial_k, initial_c_solow, savings);

time_grid = range(0, step = step_size, length = number_of_steps+1);

using Plots
p1 = plot(time_grid, [k_path_best k_path_solow], label = ["Neoclassical GM" "Solow Model"], title="Capital Path", legend = :outerbottom)
savefig(p1, "3_capital_path.png")
p2 = plot(time_grid, [c_path_best c_path_solow], label = ["Neoclassical GM" "Solow Model"], title="Consumption Path", legend = :outerbottom)
savefig(p2, "3_consumption_path.png")
#display(plot(p1, p2, layout = (1,2)))

p = plot(p1, p2, layout = (1,2))
savefig("3_both_paths.png")