# AiidaDFTK

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mfherbst.github.io/AiidaDFTK.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mfherbst.github.io/AiidaDFTK.jl/dev/)
[![Build Status](https://github.com/mfherbst/AiidaDFTK.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mfherbst/AiidaDFTK.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/mfherbst/AiidaDFTK.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mfherbst/AiidaDFTK.jl)

Julia-side implementation of the Aiida plugin of [DFTK](https://dftk.org).

Enables running a DFTK calculation from a JSON file,
such that all computed results are again available in JSON or HDF5 for programmatic parsing.
