using MaoyeLoc

if length(ARGS)==1 && only(ARGS) in ("--help","-h")
    println("用法：julia --project=. bin/plot_locations.jl prepared_plot_directory")
    println("从 prepare_location_plot 保存的 CSV/TOML 绘制事件分布、正交深度剖面和累计计数，输出 locations.png/pdf。")
    exit()
end
length(ARGS)==1 || error("Usage: julia --project=. bin/plot_locations.jl prepared_plot_directory")
using CairoMakie
println("Figures: ",plot_location_results(only(ARGS)))
