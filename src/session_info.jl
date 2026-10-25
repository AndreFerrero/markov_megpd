using Dates, LinearAlgebra

function print_session_info()
    println("="^60)
    println("Start time:       ", now())
    println("Host:             ", gethostname())
    println("Julia version:    ", VERSION)
    println("Project:          ", Base.active_project())
    println("Script:           ", PROGRAM_FILE)
    println("PID:              ", getpid())

    # Threads
    println("Julia threads:    ", Threads.nthreads(:default),
            " (interactive: ", Threads.nthreads(:interactive), ")")
    println("BLAS threads:     ", BLAS.get_num_threads(), " (", BLAS.get_config().loaded_libs[1].libname |> basename, ")")

    # What SLURM granted
    for v in ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_CPUS_PER_TASK",
              "SLURM_MEM_PER_NODE", "SLURM_JOB_PARTITION",
              "JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
        println(rpad(v * ":", 22), get(ENV, v, "(unset)"))
    end

    # What the OS actually lets this process use
    println("CPUs on node:     ", Sys.CPU_THREADS, " (node total, not your allocation)")
    status = read("/proc/self/status", String)
    m = match(r"Cpus_allowed_list:\s*(\S+)", status)
    println("CPUs allowed:     ", m === nothing ? "?" : m[1])

    # Memory: node vs. the limit actually enforced
    println("Node memory (GB): ", round(Sys.total_memory() / 2^30, digits=1))
    for f in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes")
        if isfile(f)
            s = strip(read(f, String))
            println("Cgroup mem limit: ", s == "max" ? "max (none)" : round(parse(Int, s) / 2^30, digits=1), " GB")
            break
        end
    end

    println("Working dir:      ", pwd())
    println("="^60)
    flush(stdout)
end

