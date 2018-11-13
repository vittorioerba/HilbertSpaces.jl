"""
    Module for Exact Diagonalization techniques.
    
    TODO TYPES:
    - Abstract half spin chains and integer spin chains to allow for easier model definition
    
    TODO PERFORMANCE:
    - Benchmarks
    - Use SparseMatrix for the hamiltonian matrix
    
    TODO BUGS:
    - Need to be able to split states at any dof in the system => need to check on index of splitting and throw error
"""
#==============================================================================#

module HilbertSpaces

using LinearAlgebra

include("Utilities.jl")
using .Utilities

export 
#Types
HilbertSpace,
HeisenbergModel1D, 
#Abstract methods
basis,
half_basis,
tensor_basis,
H,
N,
#Functions
H_mat,
getGS,
split_state,
singular_values,
entanglement

#==============================================================================#
"""
    Abstract type to contain the information on a finite dimensional quantum system.
    
    T is the type of the variables that represent basis states.

    Each concrete istance M<:HilbertSpace is expected to have the following methods:
    - basis(::M) :: Array{T,1}
        Returns an array containing the basis for the Hilbert space
    - basis(::M, ::Int64) :: Array{T,1}
        Returns an array containing the basis for the Hilbert space with specified dof
    - half_basis(::M) :: Array{T,1}
        Returns the basis of the system, but with half dof (assuming even dof)
    - tensor_basis(::M,::T,::T) ::T
        Tensor product of two states
    - H(::M, A::T, B::T) :: Array{Float64,2} (forse meglio Complex64?)
        Returns the energy matrix in the basis specific to the model
    - N(::M) :: Int64
        Returns the number of dof of the system
"""

abstract type HilbertSpace{T} end

"""
    Abstract function to compute the energy matrix of a generic HilbertSpace in its chosen basis.
"""

function H_mat(M::HilbertSpace{T}) where T

    base = basis(M)
    size = length(base)
    matrix = zeros(size,size)

    for (i,b1) in enumerate(base)
        for (j,b2) in enumerate(base)
            # set upper diagonal
            if(i<=j)
                matrix[j,i] = H(M,base[j],base[i])
            # set lower diagonal
            else
                matrix[j,i] = conj(matrix[i,j])
            end
        end
    end

    return matrix
end

"""
    Abstract function to compute the groundstate of the system in the specified basis
"""

function getGS(M::HilbertSpace{T}) where T
    return eigvecs(H_mat(M))[:,1]
end

"""
    Abstract function that decomposes the basis as a tensor product of two elements, one of d dof, and the other of N(M) - d dof, and computes the matrix representation of a state defined by
    psi_{l,r} = psi_{l⨂r} for l in basis(d dof) and r in basis(N(M)-d dof)
    
    Due to type stability, the current version of the function supports only d=n/2 and even n, without checks!
"""

function split_state(M::HilbertSpace{T}, state::Array{Float64,1}, d::Int64) where T
    # should check on 0<d<N(M) and give error
    
    lb = basis(M,d) # left basis
    ldof = length(lb)
    rb = basis(M,N(M)-d) #right basis
    rdof = length(rb)
    
    psi = zeros(ldof,rdof)
    
    #loop on basis
    for (r, hb1) in enumerate(rb)
        for (l, hb2) in enumerate(lb)
            tensor = tensor_basis(M,hb2,hb1) #join the 2 basis elements with tensor product
            index = findfirst(x->x==tensor, basis(M)) #find the right element in psi
            psi[l,r] = state[index]
        end
    end
    
    return psi
    
end

# function split_state(M::HilbertSpace{T}, state::Array{Float64,1}) where T
#
#     hb = basis(M,div(N(M),2))
#     dof = length(hb)
#
#     psi = zeros(dof,dof)
#
#     for (r, hb1) in enumerate(hb)
#         for (l, hb2) in enumerate(hb)
#             tensor = tensor_basis(M,hb2,hb1)
#             index = findfirst(x->x==tensor, basis(M))
#             psi[l,r] = state[index]
#         end
#     end
#
#     return psi
# end

"""
    Singular values of a state, for a sistem split in half (Schmidt numbers)
"""

function singular_values(M::HilbertSpace{T}, state::Array{Float64,1}) where T
   return svd( split_state(M, state, div(N(M),2) ),full=true ).S    
end

"""
    Entanglement entropy of a state, for a sistem split in half
"""

function entanglement(M::HilbertSpace{T}, state::Array{Float64,1}) where T   
   return -sum(map( x -> x^2 * log(x^2) , singular_values(M,state))) 
end

#==============================================================================#

#= ABSTRACT SUBTYPES =#
#TODO include("HalfSpinChain.jl")

#= CONCRETE IMPLEMENTATIONS =#
include("HeisenbergModel1D.jl")

#==============================================================================#
end




