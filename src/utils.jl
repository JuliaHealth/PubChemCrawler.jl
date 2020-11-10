"""
    rex = atomregex(chemicalsymbol)

Create a regular expression for detecting how many atoms of type `chemicalsymbol` are in a molecular formula.

# Examples

The formula for calcium bicarbonate is Ca(HCO3)2, i.e., C2CaH2O6.

```julia
julia> match(atomregex("C"), "C2CaH2O6")
RegexMatch("C2", 1="2")

julia> match(atomregex("Ca"), "C2CaH2O6")
RegexMatch("Ca", 1="")

julia> match(atomregex("H"), "C2CaH2O6")
RegexMatch("H2", 1="2")

julia> match(atomregex("O"), "C2CaH2O6")
RegexMatch("O6", 1="6")
```

Note that the regex for `"C"` does not match `"Ca"`.
"""
atomregex(str) = Regex("$str(?![a-z])(\\d*)")

function parse_formula(str)
    prs = Pair{String,Int}[]
    idx, l = 1, length(str)
    while idx <= l
        atom = str[idx:idx]
        idx = nextind(str, idx)
        if idx <= l
            c = str[idx]
            if islowercase(c)
                atom *= c
                idx = nextind(str, idx)
            end
        end
        idxn = idx
        while idx <= l && isnumeric(str[idx])
            idx = nextind(str, idx)
        end
        push!(prs, atom => idx == idxn ? 1 : parse(Int, str[idxn:idx-1]))
    end
    return prs
end

function canonicalize_properties(props)
    endswith(props, ',') && (props = props[1:end-1])
    endswith(props, '/') || (props *= '/')
    return props
end
