# Native TOML configuration and physical resources; observation files remain references.
function _state_pack(value)
    if value isa AbstractMatrix || value isa AbstractArray && ndims(value)>2
        return Dict("__maoyeloc_array__"=>true,"dimensions"=>collect(size(value)),
            "values"=>[_state_pack(x) for x in vec(value)])
    elseif value isa NamedTuple || value isa AbstractDict
        return Dict(string(k)=>_state_pack(v) for (k,v) in pairs(value))
    elseif value isa Tuple || value isa AbstractVector
        return [_state_pack(x) for x in value]
    elseif value===nothing
        return Dict("__maoyeloc_nothing__"=>true)
    elseif value===missing
        return Dict("__maoyeloc_missing__"=>true)
    elseif value isa Symbol
        return Dict("__maoyeloc_symbol__"=>String(value))
    end
    value
end

function _state_unpack(value)
    if value isa AbstractDict
        get(value,"__maoyeloc_nothing__",false) && return nothing
        get(value,"__maoyeloc_missing__",false) && return missing
        haskey(value,"__maoyeloc_symbol__") && return Symbol(value["__maoyeloc_symbol__"])
        if get(value,"__maoyeloc_array__",false)
            return reshape([_state_unpack(x) for x in value["values"]],Tuple(value["dimensions"]))
        end
        return Dict{String,Any}(k=>_state_unpack(v) for (k,v) in value)
    elseif value isa AbstractVector
        return [_state_unpack(x) for x in value]
    end
    value
end

function _state_resources!(config,base)
    resources=Dict{String,Any}[]
    function capture!(table,key,role)
        haskey(table,key) || return
        original=table[key]
        original isa AbstractString || return
        absolute=normpath(joinpath(base,original))
        isfile(absolute) || throw(ArgumentError("Physical state resource not found: $absolute"))
        name="physical_$(length(resources)+1)"*splitext(absolute)[2]
        push!(resources,Dict("name"=>name,"role"=>role,"source_path"=>absolute,"contents"=>read(absolute,String)))
        table[key]="__state_resource__/"*name
    end
    function model!(table)
        table isa AbstractDict || return
        for key in ("p_table","s_table"); capture!(table,key,"travel_time_table"); end
        for key in ("p_model","s_model","ray_model"); haskey(table,key) && model!(table[key]); end
        for child in values(get(table,"models",Dict())); model!(child); end
    end
    haskey(config,"model") && model!(config["model"])
    files=get(config,"files",Dict())
    capture!(files,"stations","stations")
    # Keep observation paths as references; a station/model state does not
    # silently absorb real phases, event catalogs or differential data.
    for (key,value) in files
        value isa AbstractString && !startswith(value,"__state_resource__/") && (files[key]=normpath(joinpath(base,value)))
    end
    value=get(config,"processing_metadata",nothing)
    value isa AbstractString && (config["processing_metadata"]=normpath(joinpath(base,value)))
    resources
end

"""
    save_physical_state(path, config_or_path; base_directory=pwd(), metadata=nothing)

Save the complete native configuration, including nested model/phase/regional
selectors, delays, weights, histories, magnitude definitions and processing
controls. Embed referenced station and native travel-time table files as text.
Observation/event files remain explicit absolute references and are not copied.
No Julia or Fortran memory image is used. The native state retains Float64
numbers and missing/empty distinctions; optional metadata retains provenance.
"""
function save_physical_state(path,config_or_path;base_directory=pwd(),metadata=nothing)
    config=config_or_path isa AbstractString ? TOML.parsefile(config_or_path) : deepcopy(config_or_path)
    config isa AbstractDict || throw(ArgumentError("A native state needs a complete configuration table or TOML path."))
    base=config_or_path isa AbstractString ? dirname(abspath(config_or_path)) : abspath(base_directory)
    resources=_state_resources!(config,base)
    state=Dict("format"=>"MaoyeLoc physical state 2","configuration"=>_state_pack(config),
        "resources"=>resources,"metadata"=>_state_pack(metadata),
        "scope"=>"Complete native physical configuration; observation files are references, not embedded data.")
    open(io->TOML.print(io,state),path,"w")
    (;path=String(path),resource_count=length(resources))
end

"""
    load_physical_state(path; output_directory=splitext(path)[1]*"_restored")

Restore the full native configuration and embedded physical resources to a
directory. Return configuration, runnable config_path and original metadata.
The external event/phase references remain the caller's original paths.
The format is MaoyeLoc physical state 2; Julia symbols, array dimensions and
missing metadata values are preserved without using a language memory image.
"""
function load_physical_state(path;output_directory=splitext(abspath(path))[1]*"_restored")
    state=TOML.parsefile(path)
    get(state,"format",nothing)=="MaoyeLoc physical state 2" ||
        throw(ArgumentError("Unknown native physical-state format."))
    configuration=_state_unpack(state["configuration"])
    output=abspath(output_directory)
    mkpath(output)
    for resource in state["resources"]
        write(joinpath(output,resource["name"]),resource["contents"])
    end
    function restore(value)
        value isa AbstractDict && return Dict(key=>restore(item) for (key,item) in value)
        value isa AbstractVector && return [restore(item) for item in value]
        value isa AbstractString && startswith(value,"__state_resource__/") && return joinpath(output,split(value,'/')[end])
        value
    end
    config=restore(configuration)
    config_path=joinpath(output,"case.toml")
    open(io->TOML.print(io,config),config_path,"w")
    (;config,config_path,output_directory=output,metadata=_state_unpack(state["metadata"]))
end

