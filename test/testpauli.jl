using LinearAlgebra  # for checking matrix properties

# Assuming MasterEquations.jl is in your current directory or installed
using MasterEquations
using MasterEquations.Pauli

# Define our test rate functions
γ = 1.0  # emission rate
w(n, m) = n*(n-1)*m     # absorption rate
v(n, m) = γ*m          # emission rate

#= # Test operator application to rate function
println("\nRate function test:")
state = [0, 1]  # test state (n=1, m=2)
ops = Eplus(1) * Eplus(1) * Eminus(2)
shifted_w, final_state = ops(w)(state...)
println("w($state) = $(w(state...))")
println("(E₁⁺E₁⁺E₂⁻w)($state) = $shifted_w at state $final_state")

# Test state-index mapping
println("\nState indexing test:")
ranges = [0:2, 0:2]
dims = length.(ranges)

# Test all state mappings
for index in 1:prod(dims)
    state = index_to_state(index, dims, ranges)
    back_index = state_to_index(state, dims, ranges)
    println("Index $index <-> State $state (back to index: $back_index)")
end =#

#= # Create transitions
ops = Eplus(1) * Eplus(1) * Eminus(2)
shifted_w = ops(w)

# Build and display rate matrix
rates = build_rate_matrix(ranges, shifted_w)
println("\nRate matrix (9x9):")
display(rates) =#

#= # Test specific transitions
println("\nSpecific transitions:")
state = [2,2] # corresponds to index 5
rate, final_state = shifted_w(state...)
i = state_to_index(final_state, dims, ranges)
j = state_to_index(state, dims, ranges)
if !isnothing(i) && !isnothing(j)
    println("From state $state (index $j) to $final_state (index $i)")
    println("Rate: $rate")
else
    println("Transition leads outside valid range")
end =#

# Create process with multiple transitions
process = step_process(
    [0:2, 0:2],
    (Eplus(1)*Eplus(1)*Eminus(2))(w),  # absorption
    (Eminus(1)*Eminus(1)*Eplus(2))(w), # inverse absorption
    (Eminus(1))(v)                      # emission
)

# Build complete rate matrix
rates = build_rate_matrix(process)#, check_conservation=true)

# Display the rate matrix
println("\nRate matrix (9x9):")
display(rates)