"""
    read_model(case_toml_path)

Read `[model]` from a case TOML file. `kind="constant_velocity_layers"` gives a
`HomogeneousModel` for one layer or a `LayeredModel` for multiple layers.
`kind="linear_gradient"` gives a `LinearGradientModel`; its scalar fields are
`top_km`, `vp_at_top_km_s`, `vs_at_top_km_s`, `gradient_p_s_inv` and
`gradient_s_s_inv`. Speeds are at the model top, not necessarily sea level.
`kind="finite_linear_gradient"` adds `bottom_km`, `vp_halfspace_km_s` and
`vs_halfspace_km_s` for a finite gradient layer over a constant half-space.
`kind="piecewise_linear_gradient"` gives continuous `depth_km`, `vp_km_s`,
`vs_km_s` node arrays and a constant half-space below the last node.
`kind="travel_time_tables"` reads a native TOML `p_table` plus explicit `vp_vs_ratio`
or `s_table`. Paths are relative to this case; `station_plane_km` is explicit.
`kind="phase_models"` supplies a `p_model` ray/table definition, plus either
positive dimensionless `s_scale` and optional `s_offset_s` in seconds
(`T_S=s_scale*T_P+s_offset_s`), or an independent `s_model` ray/table definition.
All depths use the same vertical reference as the source and station data.
`kind="station_models"` selects a fixed model for each station using
`default_model`, `[model.models.<id>]` model definitions and
`[model.station_models]` station mappings.
Optional `[model.station_delays_s.<model_id>.<station_id>]` tables hold known
`P` and/or `S` delays in seconds, added to that model's ray time.
`kind="regional_models"` adds ordered `nodes` (model_id, x_km, y_km,
inner_radius_km, outer_radius_km), default_model and models branch tables.
Each branch may itself be a fixed station_models set.
"""
function read_model(case_toml_path::AbstractString)
    config = TOML.parsefile(case_toml_path)
    model = get(config, "model", nothing)
    model isa AbstractDict || throw(ArgumentError("Missing [model] in $case_toml_path."))
    return _read_model_table(model; base=dirname(abspath(case_toml_path)))
end

function _read_table_definition(value,base)
    value isa AbstractString && return read_travel_time_table(normpath(joinpath(base,value)))
    value isa AbstractDict || throw(ArgumentError("A travel-time table is a native TOML path or table values."))
    rows=value["time_s"]
    times=reduce(vcat,(permutedims(Float64.(row)) for row in rows))
    TravelTimeTable(value["distance_km"],value["depth_km"],times;
        title=get(value,"title",""),node_error_s=get(value,"node_error_s",0.))
end

function _read_model_table(model::AbstractDict; base=pwd())
    if get(model,"kind","")=="ray_transform"
        ray=_read_model_table(model["ray_model"];base)
        return RayTransformModel(ray;source_depth_shift_km=get(model,"source_depth_shift_km",0.),
            receiver_depth_shift_km=get(model,"receiver_depth_shift_km",0.),
            receiver_plane_km=get(model,"receiver_plane_km",nothing),
            time_scale=get(model,"time_scale",1.),time_offset_s=get(model,"time_offset_s",0.))
    end
    if get(model,"kind","") == "regional_models"
        definitions = get(model,"models",nothing)
        nodes = get(model,"nodes",Any[])
        default = get(model,"default_model",nothing)
        definitions isa AbstractDict && nodes isa AbstractVector && default isa AbstractString ||
            throw(ArgumentError("Regional models require models, default_model and optional ordered nodes."))
        branches = Dict(id=>_read_model_table(table;base) for (id,table) in definitions)
        circles = RegionalNode[]
        fields = ("model_id","x_km","y_km","inner_radius_km","outer_radius_km")
        for node in nodes
            node isa AbstractDict && all(key -> haskey(node,key),fields) ||
                throw(ArgumentError("Each regional node needs "*join(fields,", ")*"."))
            push!(circles,RegionalNode((node[key] for key in fields)...))
        end
        return RegionalModels(branches,circles;default_model_id=default)
    end
    if get(model,"kind","") == "station_models"
        definitions = get(model,"models",nothing)
        assignments = get(model,"station_models",Dict{String,String}())
        default = get(model,"default_model",nothing)
        definitions isa AbstractDict && assignments isa AbstractDict && default isa AbstractString ||
            throw(ArgumentError("Station models need models, default_model and an optional station_models table."))
        leaves = Dict{String,StationModelLeaf}()
        for (id, table) in definitions
            table isa AbstractDict || throw(ArgumentError("Model $id must be a table defining a supported one-dimensional ray, travel-time lookup, or P/S relation."))
            leaves[String(id)] = _read_model_table(table; base)
        end
        delay_tables = get(model,"station_delays_s",Dict())
        delay_tables isa AbstractDict || throw(ArgumentError("station_delays_s must be model/station tables."))
        delays = Dict{Tuple{String,String,String},Float64}()
        for (model_id,stations) in delay_tables
            stations isa AbstractDict || throw(ArgumentError("Delay model $model_id must contain station tables."))
            haskey(leaves,model_id) || throw(ArgumentError("Unknown delay model ID $model_id."))
            for (station_id,phases) in stations
                phases isa AbstractDict || throw(ArgumentError("Delay station $station_id must contain P/S delays."))
                for (phase,value) in phases
                    value isa Real && isfinite(value) || throw(ArgumentError("Station delays must be finite seconds."))
                    delays[(String(model_id),String(station_id),String(phase))] = Float64(value)
                end
            end
        end
        return StationModels(leaves,assignments;default_model_id=default,station_delays_s=delays)
    end
    kind = get(model, "kind", "")
    if kind == "phase_models"
        get(model,"p_model",nothing) isa AbstractDict || throw(ArgumentError("phase_models requires a p_model table."))
        p = _read_model_table(model["p_model"];base)
        p isa PhaseModelLeaf || throw(ArgumentError("P/S components must be direct ray/table definitions, not selectors or nested phase pairs."))
        if haskey(model,"s_model")
            model["s_model"] isa AbstractDict || throw(ArgumentError("s_model must be a model table."))
            !haskey(model,"s_scale") && !haskey(model,"s_offset_s") ||
                throw(ArgumentError("An independent S model does not also use an affine S relation."))
            s = _read_model_table(model["s_model"];base)
            s isa PhaseModelLeaf || throw(ArgumentError("S component must be a direct ray/table definition."))
            return PhaseModels(p,s)
        end
        haskey(model,"s_scale") || throw(ArgumentError("An affine S relation requires explicit s_scale."))
        return PhaseModels(p;s_scale=model["s_scale"],s_offset_s=get(model,"s_offset_s",0.0))
    end
    if kind == "travel_time_tables"
        haskey(model,"p_table") || throw(ArgumentError("Travel-time tables need p_table."))
        xor(haskey(model,"s_table"),haskey(model,"vp_vs_ratio")) ||
            throw(ArgumentError("Supply exactly one of s_table and vp_vs_ratio."))
        p = _read_table_definition(model["p_table"],base)
        plane = get(model,"station_plane_km",0.0)
        return haskey(model,"s_table") ?
            TravelTimeTableModel(p,_read_table_definition(model["s_table"],base);station_plane_km=plane) :
            TravelTimeTableModel(p;vp_vs_ratio=model["vp_vs_ratio"],station_plane_km=plane)
    end
    if kind == "piecewise_linear_gradient"
        keys = ("depth_km", "vp_km_s", "vs_km_s")
        all(key -> haskey(model,key) && model[key] isa AbstractVector,keys) ||
            throw(ArgumentError("Piecewise gradients need depth_km, vp_km_s and vs_km_s arrays."))
        return PiecewiseGradientModel((model[key] for key in keys)...)
    end
    if kind == "finite_linear_gradient"
        keys = ("vp_at_top_km_s", "vs_at_top_km_s", "gradient_p_s_inv",
                "gradient_s_s_inv", "vp_halfspace_km_s", "vs_halfspace_km_s",
                "top_km", "bottom_km")
        all(key -> haskey(model,key) && model[key] isa Real, keys) ||
            throw(ArgumentError("Finite linear gradient needs scalar "*join(keys,", ")*"."))
        vp,vs,gp,gs,vph,vsh,top,bottom = (model[key] for key in keys)
        return FiniteGradientModel(vp,vs,gp,gs,vph,vsh;top_km=top,bottom_km=bottom)
    end
    if kind == "linear_gradient"
        keys = ("vp_at_top_km_s", "vs_at_top_km_s", "gradient_p_s_inv",
                "gradient_s_s_inv", "top_km")
        all(key -> haskey(model,key) && model[key] isa Real, keys) ||
            throw(ArgumentError("Linear gradient needs scalar "*join(keys,", ")*"."))
        vp, vs, gp, gs, top = (model[key] for key in keys)
        return LinearGradientModel(vp, vs, gp, gs; top_km=top)
    end
    kind == "constant_velocity_layers" ||
        throw(ArgumentError("Supported model definitions are constant_velocity_layers, linear_gradient, finite_linear_gradient, piecewise_linear_gradient, travel_time_tables and phase_models."))
    keys = ("layer_top_km", "vp_km_s", "vs_km_s")
    all(key -> haskey(model, key) && model[key] isa AbstractVector, keys) ||
        throw(ArgumentError("Model needs layer_top_km, vp_km_s and vs_km_s arrays."))
    top, vp, vs = (model[key] for key in keys)
    if length(top) == length(vp) == length(vs) == 1
        return HomogeneousModel(only(vp), only(vs); top_km=only(top))
    end
    return LayeredModel(top, vp, vs)
end

# Read one CSV record at a time. A quoted field may contain commas, escaped
# quotes and newlines; its observation identity uses the first physical line.
function _csv_record(io, path, line_number; strip_unquoted=true)
    while !eof(io)
        raw = readline(io; keep=true)
        line_number[] += 1
        isempty(strip(raw)) && continue
        first_line = line_number[]
        if !occursin('"', raw)
            # split on a String, as in the original unquoted reader. Splitting
            # a chomp-created SubString adds iterator-state allocations per field.
            fields = split(String(chomp(raw)), ','; keepempty=true)
            return first_line, strip_unquoted ? String.(strip.(fields)) : String.(fields)
        end
        fields = String[]
        field = IOBuffer()
        state = :start
        quoted = false
        while true
            text = chomp(raw)
            for char in text
                if state == :quoted
                    if char == '"'
                        state = :closed
                    else
                        write(field, char)
                    end
                elseif state == :closed && char == '"'
                    write(field, char)
                    state = :quoted
                elseif char == ','
                    value = String(take!(field))
                    push!(fields, quoted || !strip_unquoted ? value : String(strip(value)))
                    state = :start
                    quoted = false
                elseif state in (:closed, :padding)
                    char in (' ', '\t') || throw(ArgumentError("$path:$first_line: unexpected text after a quoted CSV field."))
                    state = :padding
                elseif state == :start && char == '"'
                    state = :quoted
                    quoted = true
                elseif state == :start && strip_unquoted && char in (' ', '\t')
                    continue
                else
                    char == '"' && throw(ArgumentError("$path:$first_line: a CSV quote must begin a field."))
                    write(field, char)
                    state = :unquoted
                end
            end
            if state != :quoted
                value = String(take!(field))
                push!(fields, quoted || !strip_unquoted ? value : String(strip(value)))
                return first_line, fields
            end
            # Preserve LF/CRLF inside a quoted field, as written by the exporter.
            write(field, SubString(raw, ncodeunits(text)+1))
            eof(io) && throw(ArgumentError("$path:$first_line: unterminated quoted CSV field."))
            raw = readline(io; keep=true)
            line_number[] += 1
        end
    end
    return nothing
end

# The callback executes inside open, including when a row fails validation.
# Native tables map names once; existing exchange callers can keep exact order.
function _foreach_csv_row(f, path, required; optional=(), allow_extra=false,
                          allow_empty=false, reorder=true, strip_unquoted=true,
                          on_header=nothing)
    open(path) do io
        line_number = Ref(0)
        first_record = _csv_record(io, path, line_number)
        first_record === nothing && throw(ArgumentError("$path: missing CSV header."))
        line, header = first_record
        length(unique(header)) == length(header) ||
            throw(ArgumentError("$path:$line: duplicate CSV column name."))
        all(name in header for name in required) ||
            throw(ArgumentError("$path:$line: required CSV columns are "*join(required, ", ")*"."))
        allow_extra || all(name in required || name in optional for name in header) ||
            throw(ArgumentError("$path:$line: unexpected CSV column."))
        reorder || header == collect(required) ||
            throw(ArgumentError("$path:$line: unexpected CSV header."))
        names = reorder ? [String.(required)..., [name for name in header if !(name in required)]...] : header
        order = Int[findfirst(==(name), header) for name in names]
        same_order = all(index == order[index] for index in eachindex(order))
        on_header === nothing || on_header(names)
        while true
            record = _csv_record(io, path, line_number; strip_unquoted)
            record === nothing && break
            line, fields = record
            length(fields) == length(header) ||
                throw(ArgumentError("$path:$line: expected $(length(header)) columns."))
            allow_empty || all(!isempty, fields) || throw(ArgumentError("$path:$line: empty field."))
            f(line, same_order ? fields : fields[order], names)
        end
    end
    return nothing
end

# Retain the materialized helper and exact header contract for exchange readers.
function csv_rows(path, header)
    rows = Tuple{Int,Vector{String}}[]
    _foreach_csv_row(path, header; reorder=false) do line, fields, names
        push!(rows, (line, fields))
    end
    return rows
end

function csv_number(text, path, line, field)
    value = tryparse(Float64, text)
    value !== nothing && isfinite(value) ||
        throw(ArgumentError("$path:$line: $field must be a finite number."))
    return value
end

function csv_sigma(text, path, line)
    sigma = csv_number(text, path, line, "sigma_s")
    sigma > 0 || throw(ArgumentError("$path:$line: sigma_s must be positive."))
    return sigma
end

function csv_phase(phase, path, line)
    phase in ("P", "S") || throw(ArgumentError("$path:$line: phase must be P or S."))
    return phase
end

"""Read station_id,x_km,y_km,z_km, preserving station IDs as strings."""
function read_stations(path::AbstractString)
    stations = Dict{String,NTuple{3,Float64}}()
    _foreach_csv_row(path, ("station_id", "x_km", "y_km", "z_km")) do line, fields, header
        id = fields[1]
        haskey(stations, id) && throw(ArgumentError("$path:$line: duplicate station_id $id."))
        stations[id] = ntuple(k -> csv_number(fields[k+1], path, line, "coordinate"), 3)
    end
    return stations
end

"""Read event_id,x_km,y_km,z_km,origin_s without renumbering event IDs."""
function read_events(path::AbstractString)
    events = NamedTuple{(:event_id, :source_km, :origin_s),Tuple{String,NTuple{3,Float64},Float64}}[]
    ids = Set{String}()
    _foreach_csv_row(path, ("event_id", "x_km", "y_km", "z_km", "origin_s")) do line, fields, header
        id = fields[1]
        id in ids && throw(ArgumentError("$path:$line: duplicate event_id $id."))
        push!(ids, id)
        source_km = ntuple(k -> csv_number(fields[k+1], path, line, "coordinate"), 3)
        origin_s = csv_number(fields[5], path, line, "origin_s")
        push!(events, (event_id=id, source_km=source_km, origin_s=origin_s))
    end
    return events
end

"""Read the event identity, position, origin and status from a native result table.

No status is selected here. Callers explicitly choose events for later work;
fixed references, failed solves and events without observations remain distinct.
"""
function read_location_results(path::AbstractString)
    event_type = NamedTuple{(:event_id,:source_km,:origin_s,:status),Tuple{String,NTuple{3,Float64},Float64,Symbol}}
    eligible_type = NamedTuple{(:event_id,:source_km,:origin_s,:status,:catalog_eligible),Tuple{String,NTuple{3,Float64},Float64,Symbol,Bool}}
    events = event_type[]
    seen = Set{String}()
    required = ("event_id", "x_km", "y_km", "z_km", "origin_s", "status")
    _foreach_csv_row(path, required; allow_extra=true, allow_empty=true,
                     strip_unquoted=false,
                     on_header=header -> (events = "catalog_eligible" in header ? Union{event_type,eligible_type}[] : event_type[])) do line, fields, header
        id = fields[1]
        !isempty(id) && !(id in seen) || throw(ArgumentError("$path:$line: duplicate or missing event ID."))
        push!(seen, id)
        position = ntuple(k -> csv_number(fields[k+1], path, line, "coordinate"), 3)
        origin = csv_number(fields[5], path, line, "origin_s")
        isempty(fields[6]) && throw(ArgumentError("$path:$line: missing event status."))
        event = (event_id=id, source_km=position, origin_s=origin, status=Symbol(fields[6]))
        eligible_index = findfirst(==("catalog_eligible"), header)
        if eligible_index !== nothing && !isempty(fields[eligible_index])
            fields[eligible_index] in ("true", "false") || throw(ArgumentError("$path:$line: invalid catalog_eligible."))
            event = merge(event, (catalog_eligible=fields[eligible_index]=="true",))
        end
        push!(events, event)
    end
    return events
end

"""Read absolute picks, preserving optional observation_id or assigning one per source row.
The ID is unique within each event; repeated station/phase readings remain separate.
"""
function read_picks(path::AbstractString)
    pick_type = NamedTuple{(:event_id,:station_id,:phase,:arrival_s,:sigma_s,:observation_id),Tuple{String,String,String,Float64,Float64,String}}
    quality_type = NamedTuple{(:event_id,:station_id,:phase,:arrival_s,:sigma_s,:observation_id,:quality_code),Tuple{String,String,String,Float64,Float64,String,Int}}
    picks = pick_type[]
    base_header=("event_id", "station_id", "phase", "arrival_s", "sigma_s")
    seen = Set{Tuple{String,String}}()
    _foreach_csv_row(path, base_header; optional=("quality_code", "observation_id"),
                     on_header=header -> (picks = "quality_code" in header ? quality_type[] : pick_type[])) do line, fields, header
        quality_index = findfirst(==("quality_code"), header)
        identity_index = findfirst(==("observation_id"), header)
        event_id, station_id = fields[1:2]
        phase = csv_phase(fields[3], path, line)
        observation_id=identity_index===nothing ? "line_$line" : fields[identity_index]
        key = (event_id,observation_id)
        key in seen && throw(ArgumentError("$path:$line: duplicate observation_id within event $event_id."))
        push!(seen, key)
        arrival_s = csv_number(fields[4], path, line, "arrival_s")
        sigma_s = csv_sigma(fields[5], path, line)
        pick=(;event_id,station_id,phase,arrival_s,sigma_s,observation_id)
        if quality_index!==nothing
            quality=tryparse(Int,fields[quality_index])
            quality!==nothing && 0<=quality<=9 || throw(ArgumentError("$path:$line: quality_code must be 0 through 9."))
            pick=merge(pick,(quality_code=quality,))
        end
        push!(picks,pick)
    end
    return picks
end

"""
Read complete arrival-time differences dt_s = arrival_i - arrival_j.
CT/CC source labels and event order are retained. sigma_s is supplied by the
user; correlation coefficients and legacy HypoDD weights are not converted.
"""
function read_differentials(path::AbstractString)
    differential_type = NamedTuple{(:event_i,:event_j,:station_id,:phase,:dt_s,:sigma_s,:source),Tuple{String,String,String,String,Float64,Float64,String}}
    quality_type = NamedTuple{(:event_i,:event_j,:station_id,:phase,:dt_s,:sigma_s,:source,:original_quality),Tuple{String,String,String,String,Float64,Float64,String,Float64}}
    differentials = differential_type[]
    seen = Set{Tuple{String,String,String,String,String}}()
    base_header = ("event_i", "event_j", "station_id", "phase", "dt_s", "sigma_s", "source")
    _foreach_csv_row(path, base_header; optional=("original_quality",),
                     on_header=header -> (differentials = "original_quality" in header ? quality_type[] : differential_type[])) do line, fields, header
        event_i, event_j, station_id = fields[1:3]
        event_i != event_j || throw(ArgumentError("$path:$line: a differential needs two distinct events."))
        phase = csv_phase(fields[4], path, line)
        source = fields[7]
        source in ("CT", "CC") || throw(ArgumentError("$path:$line: source must be CT or CC."))
        # Reversing a pair describes the same measurement, not new information.
        first_id, second_id = minmax(event_i, event_j)
        key = (first_id, second_id, station_id, phase, source)
        key in seen && throw(ArgumentError("$path:$line: duplicate pair/station/phase/source differential."))
        push!(seen, key)
        dt_s = csv_number(fields[5], path, line, "dt_s")
        sigma_s = csv_sigma(fields[6], path, line)
        row=(;event_i,event_j,station_id,phase,dt_s,sigma_s,source)
        length(header)==8 && (row=merge(row,(original_quality=csv_number(fields[8],path,line,"original_quality"),)))
        push!(differentials,row)
    end
    return differentials
end

# Native in-memory inputs use the numerical records accepted by locate/relocate.
# Keep extra fields and input order; only normalize the fields read from CSV.
function _memory_id(value,label)
    value isa AbstractString && !isempty(strip(value)) ||
        throw(ArgumentError("$label must be a nonempty string; numeric IDs are not inferred."))
    String(value)
end

function _memory_number(value,label)
    value isa Real && !(value isa Bool) && isfinite(value) ||
        throw(ArgumentError("$label must be a finite number."))
    number=Float64(value)
    isfinite(number) || throw(ArgumentError("$label is outside the Float64 range."))
    number
end

function _memory_xyz(value,label)
    (value isa Tuple || value isa AbstractVector) && length(value)==3 ||
        throw(ArgumentError("$label must contain x, y, z in km."))
    ntuple(index->_memory_number(value[index],label),3)
end

function _memory_rows(rows,fields,label)
    rows isa AbstractVector && all(row->row isa NamedTuple && all(key->haskey(row,key),fields),rows) ||
        throw(ArgumentError("inputs.$label must be a vector of NamedTuple records containing $(join(fields, ", "))."))
    rows
end

function _memory_stations(data)
    data isa AbstractDict || throw(ArgumentError("inputs.stations must map string station IDs to (x, y, z) in km."))
    Dict{String,NTuple{3,Float64}}(_memory_id(id,"station_id")=>_memory_xyz(xyz,"Station $id coordinates") for (id,xyz) in data)
end

function _memory_events(data)
    rows=_memory_rows(data,(:event_id,:source_km,:origin_s),"events")
    isempty(rows) && return NamedTuple[]
    seen=Set{String}()
    map(enumerate(rows)) do (index,row)
        event_id=_memory_id(row.event_id,"Memory event $index event_id")
        event_id in seen && throw(ArgumentError("Duplicate memory event_id $event_id."))
        push!(seen,event_id)
        merge(row,(;event_id,source_km=_memory_xyz(row.source_km,"Event $event_id source_km"),
            origin_s=_memory_number(row.origin_s,"Event $event_id origin_s")))
    end
end

function _memory_picks(data)
    rows=_memory_rows(data,(:event_id,:station_id,:phase,:arrival_s,:sigma_s),"picks")
    isempty(rows) && return NamedTuple[]
    seen=Set{Tuple{String,String}}()
    map(enumerate(rows)) do (index,row)
        label="Memory pick $index"
        event_id=_memory_id(row.event_id,"$label event_id")
        station_id=_memory_id(row.station_id,"$label station_id")
        phase=csv_phase(_memory_id(row.phase,"$label phase"),"memory picks",index)
        observation_id=haskey(row,:observation_id) ? _memory_id(row.observation_id,"$label observation_id") : "memory_row_$index"
        key=(event_id,observation_id)
        key in seen && throw(ArgumentError("$label: duplicate observation_id within event $event_id."))
        push!(seen,key)
        sigma_s=_memory_number(row.sigma_s,"$label sigma_s")
        sigma_s>0 || throw(ArgumentError("$label sigma_s must be positive."))
        pick=merge(row,(;event_id,station_id,phase,observation_id,sigma_s,
            arrival_s=_memory_number(row.arrival_s,"$label arrival_s")))
        if haskey(row,:quality_code)
            row.quality_code isa Integer && !(row.quality_code isa Bool) && 0<=row.quality_code<=9 ||
                throw(ArgumentError("$label quality_code must be an integer from 0 through 9."))
            pick=merge(pick,(quality_code=Int(row.quality_code),))
        end
        pick
    end
end

function _memory_differentials(data)
    rows=_memory_rows(data,(:event_i,:event_j,:station_id,:phase,:dt_s,:sigma_s,:source),"differentials")
    isempty(rows) && return NamedTuple[]
    seen=Set{Tuple{String,String,String,String,String}}()
    map(enumerate(rows)) do (index,row)
        label="Memory differential $index"
        event_i=_memory_id(row.event_i,"$label event_i")
        event_j=_memory_id(row.event_j,"$label event_j")
        event_i!=event_j || throw(ArgumentError("$label needs two distinct events."))
        station_id=_memory_id(row.station_id,"$label station_id")
        phase=csv_phase(_memory_id(row.phase,"$label phase"),"memory differentials",index)
        source=_memory_id(row.source,"$label source")
        source in ("CT","CC") || throw(ArgumentError("$label source must be CT or CC."))
        first_id,second_id=minmax(event_i,event_j)
        key=(first_id,second_id,station_id,phase,source)
        key in seen && throw(ArgumentError("$label: duplicate pair/station/phase/source differential."))
        push!(seen,key)
        sigma_s=_memory_number(row.sigma_s,"$label sigma_s")
        sigma_s>0 || throw(ArgumentError("$label sigma_s must be positive."))
        difference=merge(row,(;event_i,event_j,station_id,phase,source,sigma_s,
            dt_s=_memory_number(row.dt_s,"$label dt_s")))
        haskey(row,:original_quality) && (difference=merge(difference,
            (original_quality=_memory_number(row.original_quality,"$label original_quality"),)))
        difference
    end
end
