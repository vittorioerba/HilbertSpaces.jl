module Utilities

export
#functions
count_in,
find_only_one,
xor

#==============================================================================#
"""
    Function to count occurencies of what in the string where, counting overlapping istances of what.
"""
function count_in(where_::String, what::String)
    if(length(what)==1)
        return length(filter(x->x==Char(what[1]), where_))
    end
    numfinds = 0
    starting = 1
    while true
        location = findnext(what, where_, starting)
        location==nothing && return numfinds
        numfinds += 1
        starting = location.stop
    end
end

"""
    If what is in where_ exactly one time, then returns the starting index of its posisiton.
    Otherwise returns -1.
"""
function find_only_one(where_::String, what::String)
    if(length(what)==1)
        return length(filter(x->x==Char(what[1]), where_))
    end
    numfinds = 0
    starting = 1
    last = 0
    while true
        location = findnext(what, where_, starting)
        location == nothing && (return numfinds==1 ? last : -1)
        numfinds += 1
        starting = location.stop
        last = location.start
    end
end
#==============================================================================#

end