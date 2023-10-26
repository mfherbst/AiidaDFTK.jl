parse_kwargs_value(value; kwargs...) = value
function parse_kwargs_value(value::AbstractString; interpolations::AbstractDict)
    if startswith(value, ":")
        Symbol(value[2:end])  # Interpret string as symbol
    elseif startswith(value, "\$")
        interpolations[value[2:end]]  # Lookup an interpolation
    else
        value
    end
end
function parse_kwargs_value(value::AbstractVector; interpolations::AbstractDict)
    parse_kwargs_value.(value; interpolations)
end
function parse_kwargs_value(value::AbstractDict; interpolations::AbstractDict)
    # Lookup symbol and construct object
    if startswith(value["\$symbol"], "Smearing.")
        # TODO Horrible special-casing to get smearing to work
        symbol = getproperty(DFTK.Smearing, Symbol(value["\$symbol"][10:end]))
    else
        symbol = getproperty(DFTK, Symbol(value["\$symbol"]))
    end
    args   = parse_kwargs_value.(get(value, "\$args", ()); interpolations)
    kwargs = parse_kwargs(get(value, "\$kwargs", Dict()); interpolations)
    symbol(args...; kwargs...)
end

function parse_kwargs(kwargs::AbstractDict; interpolations=Dict{String,Any}())
    parsed = Dict{Symbol, Any}()
    for (key, value) in pairs(kwargs)
        parsed[Symbol(key)] = parse_kwargs_value(value; interpolations)
    end
    parsed
end
