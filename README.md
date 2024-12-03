# AiidaDFTK

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://epfl-matmat.github.io/AiidaDFTK.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://epfl-matmat.github.io/AiidaDFTK.jl/dev/)
[![Build Status](https://github.com/epfl-matmat/AiidaDFTK.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/epfl-matmat/AiidaDFTK.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/epfl-matmat/AiidaDFTK.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/epfl-matmat/AiidaDFTK.jl)

Julia-side implementation of the [Aiida](https://aiida.net)
plugin of [DFTK](https://dftk.org).
For more information how to setup and use DFTK with Aiida,
see the [documentation of the python side](https://github.com/aiidaplugins/aiida-dftk)
of the plugin.

This package implements a `run` function, which parses a JSON input
to trigger a DFTK-based calculation.
Results are again dumped in a JSON/HDF5-compatible form
and can be easily parsing in downstream scripts.

At the moment there are no guarantees in the input / output format
except that it works in combination with the `aiida-dftk` plugin.
If you wish to use the pipeline implemented by this package
in combination with a different downstream driver,
please feel free to contact us.
