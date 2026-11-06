using MaoyeLoc, TOML

if isempty(ARGS) || ARGS in (["--help"],["-h"])
    println("用法：julia --project=. bin/generate_travel_time_table.jl model_and_grid.toml output_table.toml")
    println("输入含原生 [model] 与 [table]；table 指定 distance_km、depth_km、station_plane_km，及可选 phase/title。")
    println("输出为单相 Float64 秒表；time_s 行为源深度、列为水平距离。生成节点不保证插值或外推精度。")
    exit(0)
end
length(ARGS) == 2 || error("Usage: julia --project=. bin/generate_travel_time_table.jl model_and_grid.toml output_table.toml")
definition = TOML.parsefile(ARGS[1])["table"]
table = generate_table(read_model(ARGS[1]);distance_km=definition["distance_km"],depth_km=definition["depth_km"],
    phase=get(definition,"phase","P"),station_plane_km=definition["station_plane_km"],title=get(definition,"title",""))
write_travel_time_table(ARGS[2],table)
println("Travel-time table: ",abspath(ARGS[2]))
