# We use JLD2 instead of HDF5 because of linking issues with libhdf5 on some systems when running through MPI:
# - HDF5 needs to be configured to use a libhdf5 implementation compatible with the used MPI.
# - On some systems, the system libhdf5 that matches the system MPI requires a newer libcurl than what Julia itself loads,
#   leading to a dynamic linking error (lib/julia/libcurl.so: version `CURL_OPENSSL_4' not found (required by /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so))
#   that can only be worked around using LD_PRELOAD.
using JLD2

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
        jldopen(filename, "w") do f
            store_hdf5_master(f, "", data)
        end
    end
end
