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
This package implements a `run` function, which parses a JSON input
to trigger a DFTK-based calculation.
Results are again dumped in a JSON/HDF5-compatible form
and can be easily parsing in downstream scripts.

At the moment there are no guarantees in the input / output format
except that it works in combination with the `aiida-dftk` plugin.
If you wish to use the pipeline implemented by this package
in combination with a different downstream driver,
please feel free to contact us.

Some general (and not necessarily up to date) remarks
on the input and output formats are given in [Input and Output format](@ref).

## API and supported methods
```@autodocs
Modules = [AiidaDFTK]
```
