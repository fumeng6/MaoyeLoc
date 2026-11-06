using MaoyeLoc
if length(ARGS)==1 && only(ARGS) in ("--help","-h")
    println("用法：julia --project=. bin/review_catalog.jl review.toml")
    println("配置须给出 case、output_directory 和 steps；按事件 ID 与 observation_id 显式编辑读数。")
    println("动作包括 recalculate、accept、skip、delete；修改后须重新定位才能接受。")
    exit()
end
length(ARGS)==1 || error("Usage: julia --project=. bin/review_catalog.jl review.toml")
result=run_catalog_review(only(ARGS))
println("Saved event review decisions and edits: ",result.output_directory)
println("Accepted catalog: ",result.accepted_catalog)
println("Working input for another review: ",result.working_config)
