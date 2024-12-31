using MasterEquations
using MasterEquations.Pauli

struct ModelParameters
    K::Float64
    S::Int64
    ModelParameters(; K, S) = new(K, S)  # Simple keyword constructor
end

function create_model(p::ModelParameters)
    # Helper function
    f(s) = p.K*s/p.S
    
    # Group transition rates
    rates = (
        up   = s -> (p.S-s)*exp(f(s)),
        down = s -> (p.S+s)*exp(-f(s)),
        #free = s -> -(p.S-s)*exp(f(s)) - (p.S+s)*exp(-f(s))
    )
    
    # Group transitions with their operators
    transitions = [
        (Eplus(1))(rates.down),
        (Eminus(1))(rates.up),
        (I())(rates.up,-1),    # Split the diagonal terms
        (I())(rates.down,-1)
    ]
    
    # Create process
    proc = step_process([-p.S:p.S], transitions...)
    
    return pauli_dynamics(process=proc)
end


p = ModelParameters(K=1.5, S=25)
model = create_model(p)

display(build_rate_matrix(model.process))