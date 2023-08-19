parse_kwargs_value(value) = value
function parse_kwargs_value(value::AbstractString)
    # Interpret string as symbol
    startswith(value, ":") ? Symbol(value[2:end]) : value
end
function parse_kwargs_value(value::AbstractVector)
    parse_kwargs_value.(value)
end

function parse_kwargs(kwargs::AbstractDict)
    parsed = Dict{Symbol, Any}()
    for (key, value) in pairs(kwargs)
        if value isa AbstractDict
            # Lookup symbol and construct object
            if startswith(value["\$symbol"], "Smearing.")
                # TODO Horrible special-casing to get smearing to work
                symbol = getproperty(DFTK.Smearing, Symbol(value["\$symbol"][10:end]))
            else
                symbol = getproperty(DFTK, Symbol(value["\$symbol"]))
            end
            args   = get(value, "\$args", ())
            kwargs = parse_kwargs(get(value, "\$kwargs", Dict()))

            parsed[Symbol(key)] = symbol(args...; kwargs...)
        else
            parsed[Symbol(key)] = parse_kwargs_value(value)
        end
    end
    parsed
end
