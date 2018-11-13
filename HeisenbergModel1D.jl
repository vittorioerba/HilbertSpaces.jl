"""
    Concrete implementation of HilbertSpace for the 1D quantum Heisenberg model with transverse magnetic field h (in x direction) and couplings (Jx,Jy,Jz) for the three n.n. interactions.
    
    The explicit expression is:
    H = SUM_<ij> { Jx Sx_i Sx_j + ...} + SUM_<i> h Sx_i 
    where the first SUM is intended over all n.n, and the second SUM is over all sites.
    S is the spin operator, i.e. 1/2 sigma if sigma are the Pauli matrices.
    
    The basis chosen is the slater basis for the single particle z-directed spin basis, and is represented as the integer defined by base2 to base10 conversion of the states written in up=1, down=0.
    
    The chain has n spins.
"""

struct HeisenbergModel1D <: HilbertSpace{Int64}
    j :: Array{Float64,1}
    h :: Float64
    n :: Int64

    function HeisenbergModel1D(j::Array{Float64,1}, h::Float64, n::Int64)
        return new(j,h,n)
    end
end

"""
    The basis is simply numbered from 0 to 2^n-1, as the bitstring (properly cut) is already the desired two spin representation
"""
function basis(M::HeisenbergModel1D) 
    return collect(0:2^M.n-1)
end

"""
    The basis is simply numbered from 0 to 2^n, assuming n even
"""
function basis(M::HeisenbergModel1D, n::Int64) 
    return collect(0:2^n-1)
end

"""
    The tensor product is the concatenation of the bit representation (without all the trailing zeroes!) for half dof states
"""
function tensor_basis(M::HeisenbergModel1D, A::Int64, B::Int64) 
    stateA = bitstring(A)[end+1-div(M.n,2):end]
    stateB = bitstring(B)[end+1-div(M.n,2):end]
    return parse(Int, stateA * stateB, base=2)
end


"""
    Implementation of the energy.
    As the x and y Pauli matrices act on 0/1 by changing it (modulo factors), Jx and Jy contribute only for states that differ by exactly two adjacent spins ()
"""
function H(M::HeisenbergModel1D, A::Int64, B::Int64)
    # convert from counting states from 1 to representing them from 000...
    #= A -=1 =#
    #= B -=1 =#

    energy = 0
    diff = bitstring(xor(A,B))[end+1-M.n:end]
    stateA = bitstring(A)[end+1-M.n:end]

    #x and y contribute only for states that have exactly two spins adjacent spins that are 01 and 10 respectively
    if length(filter(x->x=='1',diff)) == 2
    
        loc = find_only_one(diff, "11")
        if loc >= 0 
            energy += M.j[1]
            if stateA[loc : loc+1] == "10" || stateA[loc : loc+1] == "01"
                energy += M.j[2]
            else
                energy -= M.j[2]
            end
        end
    end
    
    #z contributes only for equal states
    if A==B
        aligned = count_in(stateA,"11") + count_in(stateA,"00")
        energy += (2*aligned+1-M.n)*M.j[3]
    end
    
    energy/=4 #pairs of spin have a 1/4 factor due to the definition of spin = 1/2 pauli

    # external field in x contributes only for states with exactly one spin flipped
    if length(filter(x->x=='1', diff))==1
        energy += M.h/2 #pairs of spin have a 1/2 factor due to the definition of spin = 1/2 pauli
    end

    return energy
end   

"""
"""
function N(M::HeisenbergModel1D) 
    return M.n
end