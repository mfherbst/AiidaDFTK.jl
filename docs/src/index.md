```@meta
CurrentModule = AiidaDFTK
```

# AiidaDFTK

Julia-side implementation of the Aiida plugin of [DFTK](https://dftk.org).

Enables running a DFTK calculation from a JSON file,
such that all computed results are again available in JSON or HDF5 for programmatic parsing.

!!! note "TODO"
    We should document the expected Julia setup (with regards to MPI, MKL (Blas / FFT), etc.).

## API and supported methods

```@autodocs
Modules = [AiidaDFTK]
```
