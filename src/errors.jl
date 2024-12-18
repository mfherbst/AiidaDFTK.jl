# We write a few key messages to an error file that we can parse from
# the AiiDA plugin to detect common failure modes.

ERRORFILE = "errors.log"

ERR_IMPORTS_SUCEEDED = "Imports succeeded"
ERR_VERSION_OK = "Package version matches requirements"
ERR_VERSION_MISMATCH = "Package version mismatch"
ERR_FINISHED_SUCCESSFULLY = "Finished successfully"

function report_error(msg)
    mpi_master() || return
    open(ERRORFILE, "a") do io
        println(io, msg)
    end
    println(msg) # also print to stdout
end

function fatal_error(msg)
    report_error(msg)
    error(msg)
end
