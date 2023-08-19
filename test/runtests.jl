using AiidaDFTK
using AtomsBase
using AtomsBaseTesting
using DFTK
using JSON3
using StaticArrays
using Test
using Unitful
using UnitfulAtomic

@testset "AiidaDFTK.jl" begin
    system_silicon = periodic_system(
        [Atom(:Si, [ 1.28,  1.28,  1.28]u"bohr";
              pseudopotential="hgh/lda/si-q4", magnetic_moment=0),
         Atom(:Si, [-1.28, -1.28, -1.28]u"bohr";
              pseudopotential="hgh/lda/si-q4", magnetic_moment=0)],
        [[0.0, 5.13, 5.13],
         [5.13, 0.0, 5.13],
         [5.13, 5.13, 0.0]]u"bohr"
    )

    @testset "parse_kwargs" begin
        using AiidaDFTK: parse_kwargs

        @test parse_kwargs(Dict("a" => 4)) == Dict(:a => 4)
        @test parse_kwargs(Dict("a" => ":b")) == Dict(:a => :b)
        @test parse_kwargs(Dict("a" => [":b", ":c"])) == Dict(:a => [:b, :c])

        data = Dict("smearing" => Dict("\$symbol" => "Smearing.MethfesselPaxton",
                                       "\$args" => [2]))
        ref = Dict(:smearing => DFTK.Smearing.MethfesselPaxton(2))
        @test parse_kwargs(data) == ref
    end

    @testset "build_system" begin
        data = Dict("bounding_box" => [[0.0, 5.13, 5.13],
                                       [5.13, 0.0, 5.13],
                                       [5.13, 5.13, 0.0]],
                    "atoms" => [Dict("symbol" => "Si",
                                     "position" => [1.28, 1.28, 1.28],
                                     "pseudopotential" => "hgh/lda/si-q4"),
                                Dict("symbol" => "Si",
                                     "position" => [-1.28, -1.28, -1.28],
                                     "pseudopotential" => "hgh/lda/si-q4")])
        res = AiidaDFTK.build_system(Dict("periodic_system" => data))
        test_approx_eq(res, system_silicon)
    end

    @testset "build_basis" begin
        smearing = Dict("\$symbol" => "Smearing.MethfesselPaxton", "\$args" => [1])
        data = Dict(
            "model_kwargs" => Dict("xc" => [":lda_x", ":lda_c_pw"],
                                   "temperature" => 1e-3,
                                   "smearing" => smearing),
            "basis_kwargs" => Dict("kgrid" => [4, 4, 4], "Ecut" => 15),
        )
        basis = AiidaDFTK.build_basis(data, system_silicon)

        ref_model = model_DFT(system_silicon, [:lda_x, :lda_c_pw];
                              temperature=1e-3, smearing=Smearing.MethfesselPaxton(1))
        ref_basis = PlaneWaveBasis(ref_model; kgrid=[4, 4, 4], Ecut=15)

        # TODO Add == implementation for models to DFTK
        @test basis.model.model_name        == ref_model.model_name
        @test basis.model.lattice           == ref_model.lattice
        @test basis.model.n_electrons       == ref_model.n_electrons
        @test basis.model.spin_polarization == ref_model.spin_polarization
        @test basis.model.smearing          == ref_model.smearing
        @test basis.model.temperature       == ref_model.temperature
        @test basis.model.positions         == ref_model.positions
        @test basis.model.symmetries        == ref_model.symmetries

        #TODO Add == implementation for PlaneWaveBasis to DFTK
        @test basis.Ecut           == ref_basis.Ecut
        @test basis.kgrid          == ref_basis.kgrid
        @test basis.kshift         == ref_basis.kshift
        @test basis.fft_size       == ref_basis.fft_size
        @test basis.variational    == ref_basis.variational
        @test basis.kcoords_global == ref_basis.kcoords_global
        @test basis.kweights       == ref_basis.kweights
        @test basis.symmetries_respect_rgrid == ref_basis.symmetries_respect_rgrid
    end

    @testset "store_hdf5" begin
        using AiidaDFTK: store_hdf5
        using HDF5

        mktempdir() do dir
            cd(dir) do
                data = (a=1, b=1.2, c=(d=3, e=4))
                store_hdf5("test.hdf5", data)
                h5open("test.hdf5", "r") do file
                    @test read(file, "a")   == data.a
                    @test read(file, "b")   == data.b
                    @test read(file, "c/d") == data.c.d
                    @test read(file, "c/e") == data.c.e
                end
            end
        end

        mktempdir() do dir
            cd(dir) do
                data = (a=[[1,2],[3,4]], b=SA[1 2; 3 4], c=3)
                store_hdf5("test.hdf5", data)
                h5open("test.hdf5", "r") do file
                    @test read(file, "a") == [1 3; 2 4]
                    @test read(file, "b") == data.b
                    @test read(file, "c") == data.c
                end
            end
        end

        mktempdir() do dir
            cd(dir) do
                data3 = (funcname="test", value=[1,2,3.1])
                store_hdf5(joinpath(dir, "test3.hdf5"), data3)
                h5open(joinpath(dir, "test3.hdf5"), "r") do file
                    @test read(file["funcname"]) == data3.funcname
                    @test all(read(file["value"]) .== data3.value)
                end
            end
        end
    end

    @testset "Functionality test run_json" begin
        using AiidaDFTK
        inputfile = "iron.json"
        ref_energy = -117.153287

        # We need to do this below @__DIR__ because the iron.json contains relative paths
        mktempdir(@__DIR__) do dir
            # Run SCF and check we got all expected files
            cd(dir) do
                (; output_files) = run_json(joinpath(@__DIR__, inputfile))
                for file in output_files
                    @test isfile(file)
                end
            end

            # Check SCF converged to the right spot
            let scfres = load_scfres(joinpath(dir, "scfres.jld2"))
                @test scfres.converged
                @test abs(scfres.energies.total - ref_energy) < 1e-2
            end

            # Check the self_consistent_field.json has all expected keys
            open(joinpath(dir, "self_consistent_field.json")) do io
                data = JSON3.read(io)
                # Energy terms should sum to total, so if we sum all values:
                @test abs(data["energies"]["total"] - ref_energy) < 1e-2
                @test sum(values(data["energies"])) ≈ 2data["energies"]["total"]

                # TODO Put tests for all keys here, which are read by Aiida
            end
        end
    end
end
