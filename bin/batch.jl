using MaoyeLoc

if length(ARGS)==1 && only(ARGS) in ("--help","-h")
    println("用法：julia --project=. bin/batch.jl batch.toml")
    println("按原生 TOML 的 steps 次序设置参数、运行定位、查询或保存物理配置。")
    println("文件路径相对于 batch.toml；计算单位为 km、km/s、s。")
    exit()
end
length(ARGS)==1 || error("Usage: julia --project=. bin/batch.jl batch.toml")
result=run_batch(only(ARGS))
println("Batch results: ",result.output_directory)
