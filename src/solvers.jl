using OrdinaryDiffEq

# Core solving functionality
function solve_dynamics(method::PauliLike, system, args...)
    # Implementation for Pauli master equation
end

function solve_dynamics(method::BlochRedfieldMethod, system, args...)
    if method.secular
        pauli_method = convert_to_pauli(method)
        return solve_dynamics(pauli_method, system, args...)
    end
    # Full Bloch-Redfield implementation
end
