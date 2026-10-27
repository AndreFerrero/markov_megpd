using Base.Threads

println("Julia started")
println("Number of threads: ", nthreads())
flush(stdout)

Threads.@threads for i in 1:20
    println("Thread ", threadid(), " — Iteration ", i)
    flush(stdout)
    sleep(60)
end

println("Julia finished")
flush(stdout)
