# Independent scientific tasks keep their own work arrays; only results are
# merged in input order. Progress describes completion and never enters a solve.
function _ordered_tasks(f, items; tasks::Integer=1, progress=nothing)
    tasks>=1 || throw(ArgumentError("tasks must be a positive integer."))
    count=length(items)
    results=Vector{NamedTuple}(undef,count)
    workers=min(Int(tasks),Threads.nthreads(),count)
    completed=Ref(0)
    progress_lock=ReentrantLock()
    function evaluate(index)
        result=f(items[index])
        results[index]=result
        if progress!==nothing
            lock(progress_lock) do
                completed[]+=1
                progress(completed[],count,index,result)
            end
        end
        nothing
    end
    if workers<=1
        for index in eachindex(items)
            evaluate(index)
        end
    else
        # Claim the next input when a worker finishes. Event/component costs
        # can differ greatly; a fixed stride can leave most threads idle.
        next_index=Threads.Atomic{Int}(1)
        @sync for _ in 1:workers
            Threads.@spawn begin
                while true
                    index=Threads.atomic_add!(next_index,1)
                    index>count && break
                    evaluate(index)
                end
            end
        end
    end
    results
end
