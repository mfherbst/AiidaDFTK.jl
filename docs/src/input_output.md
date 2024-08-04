# Input and Output format

This file gives a quick and dirty overview of the input and files of this package.
No guarantees are given that this data is correct or up to date.
We only make sure it works well in combination with
the [aiida-dftk](https://github.com/aiidaplugins/aiida-dftk) plugin.

If you wish to use the same IO pipeline in combination with a different package,
please feel free to get in touch with us.

## JSON input
We follow the file [iron.json](https://github.com/mfherbst/AiidaDFTK.jl/raw/master/test/iron.json) as an example.
All units are atomic units.

The dict key `periodic_system` contains the structure:
`bounding_box` the list of lattice vectors,
`atoms` the individual atoms with their cartesian positions,
pseudopotential information and extra keys,
such as the initial magnetic moments.
This block of data is parsed by `AiidaDFTK.build_system`.

The dict keys `model_kwargs` and `basis_kwargs` list the keyword arguments
to be passed to the [`DFTK.model_DFT`](https://docs.dftk.org/stable/api/#DFTK.model_DFT)
and [`DFTK.PlaneWaveBasis`](https://docs.dftk.org/stable/api/#DFTK.PlaneWaveBasis)
functions. This data passes through the `AiidaDFTK.parse_kwargs` function,
such that some special processing is available:
- Strings prefixed with `:` are interpreted as Julia symbols
- Dicts with a key `$symbol` cause a Julia object to be constructed.
  In the iron example, this is used to specify Gaussian smearing
  by constructing a `Smearing.Gaussian()` object.
- The special strings `$basis` and `$model` will be replaced by the results
  to the calls to `DFTK.PlaneWaveBasis` and `DFTK.model_DFT`, respectively.
For the `functionals` key (Exchange-correlation functional) all keys supported
by `model_DFT`,
(that is the keys of [Libxc](https://www.tddft.org/programs/libxc/functionals/),
are available.

The dict key `scf` specifies the SCF options. With the `$function` key
set to `self_consistent_field` a one-shot SCF is run.
The `$kwargs` list the keyword arguments passed to the
[`DFTK.self_consistent_field`](https://docs.dftk.org/stable/api/#DFTK.self_consistent_field)
function. These keys again pass through `AiidaDFTK.parse_kwargs`,
such that above special processing is available.
In the `iron.json` example this is used to construct a `DFTK.ScfConvergenceEnergy`
object, setting the convergence tolerance to `1e-4` in the energy
and to switch the mixing to `DFTK.KerkerDosMixing`.

Finally the dict key `postscf` lists functions to run after the SCF has been converged.
This includes
- `compute_forces_cart`: Computation of Cartesian forces
- `compute_stresses_cart`: Computation of Cartesian stresses.

## Output files
After the SCF has been run,
the [`DFTK.save_scfres`](https://docs.dftk.org/stable/api/#DFTK.save_scfres)
function is used to store the final state in the `self_consistent_field.json`
file and in the JLD2 file specified by the `scf/checkpointfile` key.
The JSON file only contains "small" results (energies, eigenvalues etc.)
while the JLD2 file contains the density (and if `scf/save_ψ` is true)
also the wavefunction. Since JLD2 is compatible with HDF5 file, the latter
file can be read using `h5py` or similar packages.
Both the JSON and JLD2 formats are used from DFTK without modification,
so see the documentation
of [`DFTK.save_scfres`](https://docs.dftk.org/stable/api/#DFTK.save_scfres)
for more details.

Additionally timings are stored in `timings.json` and results of post-SCF
steps (such as force calculations) in `$function.hdf5`. I.e. if you run
the post-SCF step `compute_forces_cart`, then the forces are stored
in `compute_forces_cart.hdf5`.
