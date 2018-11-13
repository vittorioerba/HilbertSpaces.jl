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
