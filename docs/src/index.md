```@meta
CurrentModule = AiidaDFTK
```

# AiidaDFTK

Julia-side implementation of the [Aiida](https://aiida.net)
plugin of [DFTK](https://dftk.org).
For more information how to setup and use DFTK with Aiida,
see the [documentation of the python side](https://github.com/aiidaplugins/aiida-dftk)
of the plugin.

Useful information might also be available in the DFTK-specific guides for
[Installation](https://docs.dftk.org/stable/guide/installation/)
and [Using DFTK on compute clusters](https://docs.dftk.org/stable/tricks/compute_clusters/).

## Functionality
This package effectively implements a pipeline where appropriate
JSON-input triggers a DFTK calculation,
which is in turn dumped in a JSON/HDF5 output form.
At the moment we give no guaranteed except that it works in the context
of the `aiida-dftk` plugin.
If you wish to use such an IO pipeline in combination with a different package,
please feel free to contact us.
Some general (and not necessarily up to date) remarks
on the input and output formats are given in [Input and Output format](@ref).

## API and supported methods
```@autodocs
Modules = [AiidaDFTK]
```
