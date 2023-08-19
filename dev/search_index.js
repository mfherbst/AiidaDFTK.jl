var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AiidaDFTK","category":"page"},{"location":"#AiidaDFTK","page":"Home","title":"AiidaDFTK","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AiidaDFTK.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AiidaDFTK]","category":"page"},{"location":"#AiidaDFTK.mpi_master","page":"Home","title":"AiidaDFTK.mpi_master","text":"Helper function to check whether we are on the master process\n\n\n\n\n\n","category":"function"},{"location":"#AiidaDFTK.run-Tuple{}","page":"Home","title":"AiidaDFTK.run","text":"run()\n\nRun a DFTK calculation from a json input file. The input file name is expected to be passed as the first argument when calling Julia (i.e. it should be available via ARGS. This function is expected to be called from queuing system jobscripts, for example:\n\njulia --project -e 'using AiidaDFTK; AiidaDFTK.run()' /path/to/input/file.json\n\nIt automatically dumps a logfile $(ARGS[1]).log (i.e. name of the input file with the log extension), which contains all logging from DFTK.\n\n\n\n\n\n","category":"method"},{"location":"#AiidaDFTK.run_json-Tuple{AbstractString}","page":"Home","title":"AiidaDFTK.run_json","text":"run_json(filename::AbstractString)\n\nRun a DFTK calculation from a json input file. Output is by default written to stdout and stderr.\n\n\n\n\n\n","category":"method"}]
}
