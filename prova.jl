include("HilbertSpaces.jl")
using .HilbertSpaces

for i in collect(2:2:10)
    H1 = HeisenbergModel1D([1.,1.,1.],0.,i)
    println(i, " | ", entanglement(H1,getGS(H1)))
    println("  | ", singular_values(H1,getGS(H1)))
end
