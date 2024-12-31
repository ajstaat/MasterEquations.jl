# Add StepProcess definition at the top
struct StepProcess
    ranges::Vector{UnitRange{Int}}
    transitions::Vector{Function}
end

# Simple constructor
function step_process(ranges::Vector{UnitRange{Int}}, transitions...)
    return StepProcess(ranges, collect(transitions))
end

# Method types
abstract type DynamicsMethod end
abstract type PauliLike <: DynamicsMethod end
abstract type BlochRedfieldLike <: DynamicsMethod end

struct PauliMethod <: PauliLike
    process::Union{StepProcess, Nothing}
    rates::Union{Matrix{Float64}, Nothing}
end

struct BlochRedfieldMethod <: BlochRedfieldLike
    secular::Bool
    # other parameters...
end

# Solution types
abstract type DynamicsSolution end

struct BasicSolution <: DynamicsSolution
    times::Vector{Float64}
    populations::Matrix{Float64}  # Each column is state populations at a time point
    # Could also store:
    # coherences::Matrix{ComplexF64}
    # observables::Dict{String, Vector{Float64}}
end

# Constructor functions
function pauli_dynamics(; process=nothing, rates=nothing)
    if process !== nothing
        return PauliMethod(process, nothing)
    elseif rates !== nothing
        return PauliMethod(nothing, rates)
    else
        error("Must specify either process or rates")
    end
end

function bloch_redfield_dynamics(; secular=false)
    method = BlochRedfieldMethod(secular)
    if secular
        # Convert to Pauli-like method under the hood
        return convert_to_pauli(method)
    end
    return method
end

# Utility functions
function convert_to_pauli(br::BlochRedfieldMethod)
    # Convert Bloch-Redfield to Pauli rates
    rates = compute_secular_rates(br)
    return PauliMethod(rates)
end
