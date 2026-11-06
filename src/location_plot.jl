function _native_plot_result_rows(path)
    rows=NamedTuple[]
    required=("event_id","x_km","y_km","z_km","origin_s","status")
    seen=Set{String}()
    _foreach_csv_row(path,required;allow_extra=true,allow_empty=true,strip_unquoted=false) do line,fields,header
        id=fields[1]
        !isempty(id) && !(id in seen) || throw(ArgumentError("$path:$line: duplicate or missing event ID."))
        push!(seen,id)
        position=ntuple(k->csv_number(fields[k+1],path,line,"coordinate"),3)
        origin=csv_number(fields[5],path,line,"origin_s")
        isempty(fields[6]) && throw(ArgumentError("$path:$line: missing event status."))
        optional(name)=begin
            index=findfirst(==(name),header)
            index===nothing || isempty(fields[index]) ? nothing : fields[index]
        end
        component=optional("component_id")
        eligible=optional("catalog_eligible")
        selection=optional("selection_status")
        exit_reason=optional("exit_reason")
        eligible===nothing || eligible in ("true","false") ||
            throw(ArgumentError("$path:$line: invalid catalog_eligible."))
        push!(rows,(event_id=id,source_km=position,origin_s=origin,status=Symbol(fields[6]),
            component_id=component===nothing ? nothing : parse(Int,component),
            selection_status=selection===nothing ? nothing : Symbol(selection),
            exit_reason=exit_reason===nothing ? nothing : Symbol(exit_reason),
            catalog_eligible=eligible===nothing ? nothing : eligible=="true"))
    end
    rows
end

function _native_plot_station_statistics(path::AbstractString)
    rows=NamedTuple[]
    required=("component_id","station_id","cc_count","ct_count","cc_rms_s","ct_rms_s")
    _foreach_csv_row(path,required;allow_extra=true,allow_empty=true,strip_unquoted=false) do line,fields,header
        cc,ct=parse.(Int,fields[3:4])
        cc>=0 && ct>=0 || throw(ArgumentError("$path:$line: negative station count."))
        rms(index,count)=count==0 ? nothing : csv_number(fields[index],path,line,"RMS")
        push!(rows,(component_id=isempty(fields[1]) ? nothing : parse(Int,fields[1]),station_id=fields[2],
            cc_count=cc,ct_count=ct,cc_rms_s=rms(5,cc),ct_rms_s=rms(6,ct)))
    end
    rows
end

"""
    prepare_location_plot(initial, result, stations; output_directory,
                          frame=nothing, vertical_datum=nothing, kwargs...)
    prepare_location_plot(result, stations; kwargs...)
    prepare_location_plot(initial_csv, result_csv, stations_csv; kwargs...)

Write full-precision native plot data for the Julia/CairoMakie renderer.
`initial` is an event vector, `result` is a native ABS/DD result or an event
vector, and `stations` maps complete station IDs to x/y/z km. The CSV method
reads a native input event table, native `events.csv`, and native station table.
All events remain present, including failed, unselected and exited events;
their exact status fields are retained, without choosing a scientific subset.

Native x/y/z are used directly. A supplied `LocalFrame` adds geographic
coordinates and supplies its unchanged vertical datum; no projection roundtrip
or depth conversion is applied. With no frame, x east/y north/z down stay in
the original Cartesian coordinates. Name their datum with `vertical_datum`;
if omitted, the plot explicitly reports it as not supplied.

Optional `station_statistics=".../station_statistics.csv"` retains supplied
statistics in the prepared data; default figures do not display RMS or compute
station statistics. Section azimuth is clockwise from north: for offsets
(dx,dy) from the declared centre, along=dx*sin(azimuth)+dy*cos(azimuth) and
cross=dx*cos(azimuth)-dy*sin(azimuth). Section width selects |cross|<=width/2
in the exported data only; it does not filter events in the default figure.
"""
function prepare_location_plot(initial,result,stations::AbstractDict;output_directory,
        frame=nothing,vertical_datum=nothing,station_statistics=nothing,
        section_azimuth_deg=0.,section_center_km=(0.,0.),section_width_km=nothing)
    frame===nothing || frame isa LocalFrame || throw(ArgumentError("frame must be a LocalFrame or nothing."))
    frame===nothing || vertical_datum===nothing || vertical_datum==frame.vertical_datum ||
        throw(ArgumentError("vertical_datum must match the supplied frame; plotting does not convert depths."))
    length(section_center_km)==2 && all(isfinite,section_center_km) && isfinite(section_azimuth_deg) ||
        throw(ArgumentError("Require a finite section centre and azimuth."))
    section_width_km===nothing || (isfinite(section_width_km) && section_width_km>0) ||
        throw(ArgumentError("Section width must be positive km."))
    final=hasproperty(result,:events) ? result.events : result
    isempty(initial) && isempty(final) && throw(ArgumentError("Both location catalogs are empty."))
    theta=mod(Float64(section_azimuth_deg),360.)
    function project(rows,default_status)
        projected=NamedTuple[]
        seen=Set{String}()
        for row in rows
            id=String(row.event_id)
            id in seen && throw(ArgumentError("Duplicate native event ID $id."))
            push!(seen,id)
            x,y,z=Float64.(row.source_km)
            all(isfinite,(x,y,z)) || throw(ArgumentError("Nonfinite native hypocentre for $id."))
            geographic=frame===nothing ? nothing : local_to_geographic(frame,(x,y,z))
            dx,dy=x-section_center_km[1],y-section_center_km[2]
            along,cross=dx*sind(theta)+dy*cosd(theta),dx*cosd(theta)-dy*sind(theta)
            push!(projected,(event_id=id,component_id=get(row,:component_id,nothing),
                latitude_deg=geographic===nothing ? nothing : geographic.latitude_deg,
                longitude_deg=geographic===nothing ? nothing : geographic.longitude_deg,
                x_km=x,y_km=y,depth_km=z,along_km=along,cross_km=cross,
                in_section=section_width_km===nothing || abs(cross)<=section_width_km/2,
                status=get(row,:status,default_status),selection_status=get(row,:selection_status,nothing),
                exit_reason=get(row,:exit_reason,nothing),catalog_eligible=get(row,:catalog_eligible,nothing),
                origin_s=get(row,:origin_s,nothing)))
        end
        projected
    end
    initial_rows=project(initial,:input_catalog)
    final_rows=project(final,:supplied_catalog)
    statistics=station_statistics isa AbstractString ? _native_plot_station_statistics(station_statistics) : station_statistics
    statistics===nothing && (statistics=NamedTuple[])
    station_rows=NamedTuple[]
    seen=Set{Tuple{Union{Nothing,Int},String}}()
    function add_station(row)
        id=String(row.station_id)
        key=(get(row,:component_id,nothing),id)
        key in seen && throw(ArgumentError("Duplicate native station/component $key."))
        push!(seen,key)
        haskey(stations,id) || throw(ArgumentError("Station statistics refer to unknown station $id."))
        point=stations[id]
        geographic=frame===nothing ? nothing : local_to_geographic(frame,point)
        cc,ct=get(row,:cc_count,0),get(row,:ct_count,0)
        cc_rms,ct_rms=get(row,:cc_rms_s,nothing),get(row,:ct_rms_s,nothing)
        cc>=0 && ct>=0 && (cc==0 || (cc_rms isa Real && isfinite(cc_rms) && cc_rms>=0)) &&
            (ct==0 || (ct_rms isa Real && isfinite(ct_rms) && ct_rms>=0)) ||
            throw(ArgumentError("Invalid native station counts or RMS for $key."))
        push!(station_rows,(station_id=id,component_id=key[1],
            latitude_deg=geographic===nothing ? nothing : geographic.latitude_deg,
            longitude_deg=geographic===nothing ? nothing : geographic.longitude_deg,
            x_km=point[1],y_km=point[2],cc_count=cc,ct_count=ct,
            cc_rms_s=cc==0 ? nothing : cc_rms,ct_rms_s=ct==0 ? nothing : ct_rms))
    end
    foreach(add_station,statistics)
    reported=Set(row.station_id for row in station_rows)
    for id in sort!(String.(collect(keys(stations))))
        id in reported || add_station((station_id=id,))
    end
    directory=abspath(output_directory)
    mkpath(directory)
    location_fields=(:event_id,:component_id,:latitude_deg,:longitude_deg,:x_km,:y_km,:depth_km,
        :along_km,:cross_km,:in_section,:status,:selection_status,:exit_reason,:catalog_eligible,:origin_s)
    for (name,rows) in (("initial.csv",initial_rows),("relocated.csv",final_rows))
        _write_result_csv(joinpath(directory,name),String.(location_fields),
            (Tuple(getproperty(row,key) for key in location_fields) for row in rows))
    end
    station_fields=(:station_id,:component_id,:latitude_deg,:longitude_deg,:x_km,:y_km,:cc_count,:ct_count,:cc_rms_s,:ct_rms_s)
    _write_result_csv(joinpath(directory,"stations.csv"),String.(station_fields),
        (Tuple(getproperty(row,key) for key in station_fields) for row in station_rows))
    datum=frame!==nothing ? frame.vertical_datum : vertical_datum===nothing ?
        "not supplied; unchanged native calculation z" : String(vertical_datum)
    frame_settings=Dict{String,Any}("kind"=>frame===nothing ? "native_cartesian" : "wgs84_local",
        "vertical_datum"=>datum,"source"=>"native_input_frame",
        "projection"=>frame===nothing ? "native Cartesian; km east/north" : "WGS84 local orthographic; km east/north")
    if frame!==nothing
        frame_settings["latitude_deg"]=frame.latitude_deg
        frame_settings["longitude_deg"]=frame.longitude_deg
    end
    initial_ids,final_ids=Set(row.event_id for row in initial_rows),Set(row.event_id for row in final_rows)
    settings=Dict("input"=>Dict("format"=>"native_full_precision","station_statistics_supplied"=>!isempty(statistics)),
        "frame"=>frame_settings,
        "section"=>Dict("azimuth_deg"=>theta,"center_x_km"=>Float64(section_center_km[1]),
            "center_y_km"=>Float64(section_center_km[2]),"width_km"=>section_width_km===nothing ? "all" : Float64(section_width_km)),
        "catalog"=>Dict("initial_count"=>length(initial_rows),"relocated_count"=>length(final_rows),
            "initial_only_ids"=>sort!(collect(setdiff(initial_ids,final_ids))),
            "relocated_only_ids"=>sort!(collect(setdiff(final_ids,initial_ids)))),
        "definitions"=>Dict("depth"=>"unchanged native calculation z km, positive down; no datum conversion",
            "rms"=>"unweighted sqrt(mean(residual_s^2)); CC and CT separate; missing class empty, true zero retained",
            "station_identity"=>"component_id and station_id; empty component means no membership supplied",
            "selection"=>"all supplied event rows retained; exact status fields saved; width affects exported oblique-section membership only"))
    open(io->TOML.print(io,settings),joinpath(directory,"plot.toml"),"w")
    (;output_directory=directory,initial=initial_rows,relocated=final_rows,stations=station_rows,frame,settings)
end

function prepare_location_plot(initial_csv::AbstractString,result_csv::AbstractString,stations_csv::AbstractString;kwargs...)
    prepare_location_plot(read_events(initial_csv),_native_plot_result_rows(result_csv),read_stations(stations_csv);kwargs...)
end

prepare_location_plot(result,stations::AbstractDict;kwargs...) =
    prepare_location_plot(NamedTuple[],result,stations;kwargs...)
prepare_location_plot(result_csv::AbstractString,stations_csv::AbstractString;kwargs...) =
    prepare_location_plot(NamedTuple[],_native_plot_result_rows(result_csv),read_stations(stations_csv);kwargs...)

_render_location_plot(directory) = throw(ArgumentError(
    "绘图需要 CairoMakie。请在当前 Julia 环境中安装 CairoMakie，并先执行 using CairoMakie；仅准备数据可用 prepare_location_plot。"))

"""
    plot_location_results(initial, result, stations; kwargs...)
    plot_location_results(result, stations; kwargs...)
    plot_location_results(prepared_directory)

Render `locations.png/pdf`: result-event horizontal distribution with stations,
two orthogonal depth sections, a narrow depth colorbar and cumulative count
against time since the first supplied origin. Geographic axes require an
explicit frame. Event status symbols and the calculation depth are preserved.
Install and load CairoMakie in the current Julia environment first.
The data-input methods return the same data as `prepare_location_plot`;
the directory method redraws its saved CSV/TOML data and returns the directory.
"""
function plot_location_results(initial,result,stations;kwargs...)
    prepared=prepare_location_plot(initial,result,stations;kwargs...)
    _render_location_plot(prepared.output_directory)
    prepared
end

plot_location_results(result,stations;kwargs...) = plot_location_results(NamedTuple[],result,stations;kwargs...)
plot_location_results(result_csv::AbstractString,stations_csv::AbstractString;kwargs...) =
    plot_location_results(_native_plot_result_rows(result_csv),read_stations(stations_csv);kwargs...)

plot_location_results(directory::AbstractString) = _render_location_plot(abspath(directory))
