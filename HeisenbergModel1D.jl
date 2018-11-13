"""
    Concrete implementation of HilbertSpace for the 1D quantum Heisenberg model with transverse magnetic field h (in x direction) and couplings (Jx,Jy,Jz) for the three n.n. interactions.
    The basis chosen is the slater basis for the single particle z-directed spin basis.
"""

struct HeisenbergModel1D <: HilbertSpace
    j :: Array{Float64,1}
    h :: Float64
    n :: Int64

    function HeisenbergModel1D(j::Array{Float64,1}, h::Float64, n::Int64)
        return new(j,h,n)
    end
end

function basis(M::HeisenbergModel1D) 
    return collect(0:2^M.n-1)
end

function H(M::HeisenbergModel1D, A::Int64, B::Int64)
    # convert from counting states from 1 to representing them from 000...
    #= A -=1 =#
    #= B -=1 =#

    energy = 0
    diff = bitstring(xor(A,B))

    #x and y contribute only for states that have exactly two spins adjacent spins that are 01 and 10 respectively
    if count_in(diff, "11") == 1
        # x and y contribution
        energy += (M.j[1] + M.j[2])
    end

    #z contributes only for equal states
    if A==B
        stateA = bitstring(A)[end+1-M.n:end]
        aligned = count_in(stateA,"11") + count_in(stateA,"00")
        energy += (2*aligned+1-M.n)*M.j[3]
    end

    # external field in x contributes only for states with exactly one spin flipped
    if length(filter(x->x=='1', diff))==1
        energy += M.h
    end

    return energy
end    


