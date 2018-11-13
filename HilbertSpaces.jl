module HilbertSpaces

export HilbertSpace, H_mat
#==============================================================================#
"""
    Abstract type to contain the information on a finite dimensional quantum system.

    Each concrete istance M<:HilbertSpace is expected to have the following methods:
    - basis(::M) :: Array{_,1}
        Returns an array containing the basis for the Hilbert space
    - H(::M) :: Array{Float64,2}
        Returns the energy matrix in the basis specific to the model
    - N() :: Int64
        Returns the number of dof of the system
"""

abstract type HilbertSpace end

"""
    Abstract function to compute the energy matrix of a generic HilbertSpace in its chosen basis.
"""

function H_mat(M::HilbertSpace)

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

#==============================================================================#
end