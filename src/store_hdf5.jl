using HDF5

convert_hdf5(data) = data
convert_hdf5(data::AbstractArray)   = Array(data)
convert_hdf5(data::Array{<:Number}) = data
convert_hdf5(data::AbstractVector{<:AbstractVecOrMat}) = reduce(hcat, data)

function store_hdf5_master(f, key::AbstractString, data)
    f[key] = convert_hdf5(data)
end
function store_hdf5_master(f, key::AbstractString, data::NamedTuple)
    for p in propertynames(data)
        store_hdf5_master(f, key * "/" * string(p), getproperty(data, p))
    end
end
function store_hdf5(filename::AbstractString, data::NamedTuple)
    if mpi_master()
        h5open(filename, "w") do f
            store_hdf5_master(f, "", data)
        end
    end
end
