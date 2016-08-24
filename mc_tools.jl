#=
Tools for working with Markov Chains

@author : Spencer Lyon, Zac Cranko

@date: 07/10/2014

References
----------

http://quant-econ.net/jl/finite_markov.html
=#

"""
Finite-state discrete-time Markov chain.

Methods are available that provide useful information such as the stationary
distributions, and communication and recurrent classes, and allow simulation of
state transitions.

##### Fields

- `p::AbstractMatrix` : The transition matrix. Must be square, all elements
must be nonnegative, and all rows must sum to unity.
- `state_values::AbstractVector` : Vector containing the values associated with
the states.
"""
type MarkovChain{T, TM<:AbstractMatrix, TV<:AbstractVector}
    p::TM # valid stochastic matrix
    state_values::TV

    function MarkovChain(p::AbstractMatrix, state_values)
        n, m = size(p)

        eltype(p) != T &&
            throw(ArgumentError("Types must be consistent with given matrix"))

        n != m &&
            throw(DimensionMismatch("stochastic matrix must be square"))

        minimum(p) <0 &&
            throw(ArgumentError("stochastic matrix must have nonnegative elements"))

        length(state_values) != n &&
            throw(DimensionMismatch("state_values should have $n elements"))

        return new(p, state_values)
    end
end

# Provide constructor that infers T from eltype of matrix
MarkovChain(p::AbstractMatrix, state_values=1:size(p, 1)) =
    MarkovChain{eltype(p), typeof(p), typeof(state_values)}(p, state_values)

"Number of states in the Markov chain `mc`"
n_states(mc::MarkovChain) = size(mc.p, 1)
