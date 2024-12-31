module Pauli

using ..MasterEquations: DynamicsMethod, PauliLike, PauliMethod, StepProcess, pauli_dynamics
import Base: *

struct StepOperator
    ops::Vector{Tuple{Int,Int}}  # Vector of (dimension, step) pairs
    
    # Inner constructor for vector of operators
    function StepOperator(ops::Vector{Tuple{Int,Int}})
        new(ops)
    end
    
    # Multiplication
    function Base.:*(a::StepOperator, b::StepOperator)
        new([a.ops..., b.ops...])
    end
    
    # Action on functions - now returns both value and shifted state
    function (op::StepOperator)(f::Function, c::Number=1)
        return function(state...)
            shifted_state = copy(collect(state))
            for (dim, step) in op.ops
                shifted_state[dim] += step
            end
            return c * f(shifted_state...), shifted_state
        end
    end
    
end

# Outer constructors - more convenient for users
Eplus(dim::Int) = StepOperator([(dim, 1)])
Eminus(dim::Int) = StepOperator([(dim, -1)])
I(dim::Int=1) = StepOperator(Vector{Tuple{Int,Int}}())  # Specify empty vector type

# Convert from linear index to state coordinates
function index_to_state(index::Int, dims::Vector{Int}, ranges::Vector{UnitRange{Int}})
    state = zeros(Int, length(dims))
    remaining = index - 1
    
    for i in length(dims):-1:1
        div, remaining = divrem(remaining, prod(dims[1:i-1]))
        state[i] = ranges[i][div + 1]
    end
    
    return state
end

# Convert from state coordinates to linear index
function state_to_index(state::Vector{Int}, dims::Vector{Int}, ranges::Vector{UnitRange{Int}})
    # Check if state is within bounds
    if !all(s in r for (s,r) in zip(state, ranges))
        return nothing
    end
    
    index = 1
    for i in 1:length(dims)
        offset = findfirst(x -> x == state[i], ranges[i]) - 1
        index += offset * prod(dims[1:i-1])
    end
    return index
end

# Build rate matrix using flattened indices
function build_rate_matrix(ranges::Vector{UnitRange{Int}}, transitions...)
    dims = length.(ranges)
    n = prod(dims)
    rates = zeros(n, n)
    
    for j in 1:n
        state = index_to_state(j, dims, ranges)
        
        for rate_func in transitions
            rate, final_state = rate_func(state...)
            
            # Get destination index if state is valid
            i = state_to_index(final_state, dims, ranges)
            if !isnothing(i)
                rates[i,j] += rate
                rates[j,j] -= rate
            end
        end
    end
    
    return rates
end

# Build rate matrix with optional conservation check
function build_rate_matrix(process::StepProcess; check_conservation::Bool=false)
    dims = length.(process.ranges)
    n = prod(dims)
    rates = zeros(n, n)
    
    # Build the rate matrix
    for j in 1:n
        state = index_to_state(j, dims, process.ranges)
        
        for rate_func in process.transitions
            rate, final_state = rate_func(state...)
            
            # Get destination index if state is valid
            i = state_to_index(final_state, dims, process.ranges)
            if !isnothing(i)
                rates[i,j] += rate
            end
        end
    end
    
    # Optional conservation check
    if check_conservation
        system_closed = true
        for i in 1:n
            row_sum = sum(rates[i,:])
            if !isapprox(row_sum, 0, atol=1e-10)
                state = index_to_state(i, dims, process.ranges)
                @warn "Process does not conserve probability at state $state" row_sum
                system_closed = false
            end
        end
        if system_closed
            @info "System is closed: probability is conserved at all states"
        end
    end
    
    return rates
end

export StepOperator, Eplus, Eminus, I, build_rate_matrix
export index_to_state, state_to_index
export StepProcess, step_process

end
