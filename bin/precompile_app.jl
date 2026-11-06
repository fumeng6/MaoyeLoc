# Build-only compilation workload. All inputs are temporary in-memory records;
# no research observations, example files or figures enter the executable bundle.
using MaoyeLoc, TOML

MaoyeLoc.cli_main(["--help"])
MaoyeLoc.cli_main(["--version"])
geographic_region(38., -122., 5.; net=2)

mktempdir() do workspace
    MaoyeLoc.cli_main(["--template", joinpath(workspace, "blank.toml")])
    models = [
        Dict("kind"=>"constant_velocity_layers", "layer_top_km"=>[-1.],
             "vp_km_s"=>[6.], "vs_km_s"=>[3.5]),
        Dict("kind"=>"constant_velocity_layers", "layer_top_km"=>[-1.,3.,8.],
             "vp_km_s"=>[4.8,6.,7.], "vs_km_s"=>[2.7,3.5,4.]),
        Dict("kind"=>"linear_gradient", "top_km"=>-1.,
             "vp_at_top_km_s"=>4.8, "vs_at_top_km_s"=>2.7,
             "gradient_p_s_inv"=>.12, "gradient_s_s_inv"=>.07),
        Dict("kind"=>"finite_linear_gradient", "top_km"=>-1., "bottom_km"=>8.,
             "vp_at_top_km_s"=>4.8, "vs_at_top_km_s"=>2.7,
             "gradient_p_s_inv"=>.12, "gradient_s_s_inv"=>.07,
             "vp_halfspace_km_s"=>7., "vs_halfspace_km_s"=>4.),
        Dict("kind"=>"piecewise_linear_gradient", "depth_km"=>[-1.,3.,8.],
             "vp_km_s"=>[4.8,6.,7.], "vs_km_s"=>[2.7,3.5,4.])]
    stations = Dict("S$k"=>(12*cosd(45k),12*sind(45k),0.) for k in 1:8)
    events = [(event_id="E$i", source_km=(0.1i,0.05i,5.0+0.1i), origin_s=10.0i) for i in 1:3]
    for (index, definition) in enumerate(models)
        model = MaoyeLoc._read_model_table(definition)
        picks = [(; event_id=e.event_id, station_id=sid, phase,
            arrival_s=e.origin_s+travel_time(model,e.source_km,stations[sid],phase).time_s,
            sigma_s=.01, quality_code=0, observation_id="$(sid)_$(phase)")
            for e in events for sid in sort!(collect(keys(stations))) for phase in ("P","S")]
        differentials = synthetic_differences(model,stations,events; seed=0).observations
        inputs = (; stations, events, picks, differentials)
        config = Dict{String,Any}("task"=>"absolute_then_relative",
            "time_reference_utc"=>"2026-01-01T00:00:00Z",
            "coordinates"=>"Local Cartesian: x east, y north, z down; km",
            "vertical_reference"=>"Common z=0 plane", "model"=>definition,
            "output_directory"=>joinpath(workspace,string(index)),
            "solver"=>Dict("max_iterations"=>50),
            "relative"=>Dict{String,Any}("use_absolute_picks"=>true,"damp"=>.5),
            "output"=>Dict("detail"=>"core"))
        for algorithm in ("dense","lsqr")
            config["relative"]["linear_solver"] = algorithm
            prepared = prepare_case(config; base_directory=workspace, inputs)
            check_case(prepared)
            effective_config(prepared)
            run_case(prepared; progress=MaoyeLoc._cli_progress_reporter(false))
        end
        config["output"]=Dict("diagnostics"=>true,"detail"=>"full")
        config["uncertainty"] = Dict("assume_independent_picks"=>true,
            "assume_independent_relative_observations"=>true,
            "relative_pairs"=>[["E1","E2"]],
            "formal"=>Dict("rderr_s"=>.15,"ercof"=>1.,"noise_model"=>"weighted_rows"))
        config["diagnostics"] = Dict("absolute"=>true,"quality_covariance_source"=>"formal")
        prepared = prepare_case(config; base_directory=workspace, inputs)
        run_case(prepared; progress=MaoyeLoc._cli_progress_reporter(true))
    end
end
