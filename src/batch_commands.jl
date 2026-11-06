# Ordered native calculations share one explicit TOML configuration.
function _batch_merge!(target,patch)
    for (key,value) in pairs(patch)
        if value isa AbstractDict
            current=get(target,key,nothing)
            table=current isa AbstractDict ? Dict{String,Any}(String(k)=>deepcopy(v) for (k,v) in pairs(current)) : Dict{String,Any}()
            target[key]=_batch_merge!(table,value)
        else
            target[key]=deepcopy(value)
        end
    end
    target
end

"""
    calculation_defaults()

Native iteration/step defaults; method selection is resolved by prepare_case,
exactly as for a single CLI case. A case must supply its own physical model,
observations, geometry, units, common time reference and observation sigmas.
"""
function calculation_defaults()
    Dict{String,Any}("solver"=>Dict{String,Any}("max_iterations"=>50,
        "location_tolerance_km"=>1e-6,"time_tolerance_s"=>1e-7))
end

function _batch_paths!(config,base)
    absolute(value)=normpath(joinpath(base,value))
    for key in ("output_directory","processing_metadata")
        get(config,key,nothing) isa AbstractString && (config[key]=absolute(config[key]))
    end
    for (key,value) in get(config,"files",Dict())
        value isa AbstractString && (config["files"][key]=absolute(value))
    end
    region=get(config,"geographic_regions",Dict())
    get(region,"reference_directory",nothing) isa AbstractString && (region["reference_directory"]=absolute(region["reference_directory"]))
    function model!(model)
        for key in ("p_table","s_table")
            get(model,key,nothing) isa AbstractString && (model[key]=absolute(model[key]))
        end
        for key in ("p_model","s_model","ray_model"); haskey(model,key) && model!(model[key]); end
        foreach(model!,values(get(model,"models",Dict())))
    end
    haskey(config,"model") && model!(config["model"])
    config
end

function _batch_query(config,paths)
    result=Dict{String,Any}()
    for path in paths
        keys=split(path,'.')
        value=config
        for key in keys
            value isa AbstractDict && haskey(value,key) || throw(ArgumentError("No current setting $path."))
            value=value[key]
        end
        result[String(path)]=deepcopy(value)
    end
    result
end

"""
    run_batch(path; output_directory=nothing, io=stdout)

Execute an ordered TOML `steps` list with persistent native configuration.
Top-level `case` or `configuration` supplies input. Actions are set (nested
`values`, optional dotted `remove`), case, query, run, review, command,
save_state, load_state, include, repeat, message, return and stop.

Includes resolve their own file paths and share current settings; end of file
returns to the parent. Repeat executes its ordered steps with persistent
changes. External command uses an explicit `argv` array and optional directory.
Each run saves its effective configuration and native result summary. Query
saves current settings without performing a calculation. CSV outputs retain
their normal per-run replacement behavior; use separate output directories
for calculations whose numerical products must be kept.
"""
function run_batch(path;output_directory=nothing,io=stdout)
    file=abspath(path)
    document=TOML.parsefile(file)
    _case_keys(document,("case","configuration","output_directory","steps"),"batch")
    base=dirname(file)
    output=output_directory===nothing ? normpath(joinpath(base,get(document,"output_directory","batch_results"))) : abspath(output_directory)
    mkpath(output)
    config=Dict{String,Any}()
    _batch_merge!(config,calculation_defaults())
    if haskey(document,"case")
        case=normpath(joinpath(base,document["case"]))
        _batch_merge!(config,_batch_paths!(TOML.parsefile(case),dirname(case)))
    end
    _batch_merge!(config,_batch_paths!(deepcopy(get(document,"configuration",Dict{String,Any}())),base))
    history=NamedTuple[]
    runs=NamedTuple[]
    queries=NamedTuple[]
    stopped=Ref(false)
    active_files=String[file]
    function save_config(index)
        directory=joinpath(output,"step_"*lpad(string(index),4,'0'))
        mkpath(directory)
        target=joinpath(directory,"case.toml")
        open(stream->TOML.print(stream,config),target,"w")
        target
    end
    function execute(steps,current_file,repeat_index=0)
        steps isa AbstractVector || throw(ArgumentError("Batch steps must be an ordered array of tables."))
        for step in steps
            stopped[] && break
            action=get(step,"action",nothing)
            action isa AbstractString || throw(ArgumentError("Every batch step requires action."))
            index=length(history)+1
            current_base=dirname(current_file)
            resolve(value)=normpath(joinpath(current_base,value))
            push!(history,(;index,action,file=current_file,repeat_index))
            if action=="set"
                patch=_batch_paths!(deepcopy(get(step,"values",Dict{String,Any}())),current_base)
                _batch_merge!(config,patch)
                for path in get(step,"remove",String[])
                    keys=split(path,'.'); table=config
                    for key in keys[1:end-1]
                        table isa AbstractDict && haskey(table,key) || throw(ArgumentError("Cannot remove unknown setting $path."))
                        table=table[key]
                    end
                    pop!(table,last(keys),nothing)
                end
            elseif action=="case"
                case=resolve(step["path"])
                empty!(config)
                _batch_merge!(config,calculation_defaults())
                _batch_merge!(config,_batch_paths!(TOML.parsefile(case),dirname(case)))
            elseif action=="query"
                effective=_batch_merge!(calculation_defaults(),config)
                selected=haskey(step,"keys") ? _batch_query(effective,step["keys"]) : deepcopy(effective)
                query=(;step_index=index,current=deepcopy(config),effective=selected)
                push!(queries,query)
                target=haskey(step,"path") ? resolve(step["path"]) : joinpath(output,"query_"*string(index)*".toml")
                mkpath(dirname(target)); open(stream->TOML.print(stream,_report_metadata(query)),target,"w")
                println(io,"Settings: ",target)
            elseif action=="run"
                haskey(step,"output_directory") && (config["output_directory"]=resolve(step["output_directory"]))
                actual=run_case(save_config(index))
                result=actual.absolute===nothing ? actual.relative : actual.absolute
                events=[(;event_id=event.event_id,status=event.status,source_km=event.source_km,origin_s=event.origin_s) for event in result.events]
                push!(runs,(;step_index=index,action,config_path=joinpath(output,"step_"*lpad(string(index),4,'0'),"case.toml"),
                    output_directory=actual.output_directory,events))
                for event in events; println(io,event.event_id," ",event.status); end
            elseif action=="review"
                run_catalog_review(resolve(step["path"]))
            elseif action=="command"
                argv=step["argv"]
                argv isa AbstractVector && !isempty(argv) && all(arg->arg isa AbstractString,argv) ||
                    throw(ArgumentError("External command requires a nonempty argv string array."))
                command=Cmd(Cmd(String.(argv));dir=haskey(step,"directory") ? resolve(step["directory"]) : current_base,windows_hide=true)
                print(io,read(command,String))
            elseif action=="save_state"
                target=resolve(step["path"]); mkpath(dirname(target))
                save_physical_state(target,config;base_directory=current_base,metadata=(;batch=file,step_index=index))
            elseif action=="load_state"
                restored=load_physical_state(resolve(step["path"]);output_directory=joinpath(output,"state_"*string(index)))
                empty!(config); _batch_merge!(config,_batch_paths!(restored.config,dirname(restored.config_path)))
            elseif action=="include"
                target=resolve(step["path"])
                target in active_files && throw(ArgumentError("Recursive batch include: $target."))
                nested=TOML.parsefile(target)
                push!(active_files,target)
                execute(nested["steps"],target,repeat_index)
                pop!(active_files)
            elseif action=="repeat"
                count=get(step,"count",nothing)
                count isa Integer && !(count isa Bool) && count>=0 || throw(ArgumentError("Repeat count must be a nonnegative integer."))
                for iteration in 1:count
                    stopped[] && break
                    execute(step["steps"],current_file,iteration)
                end
            elseif action=="message"
                println(io,step["text"])
            elseif action=="return"
                return
            elseif action=="stop"
                stopped[]=true
            else
                throw(ArgumentError("Unknown batch action $action."))
            end
        end
    end
    execute(document["steps"],file)
    final_path=joinpath(output,"final_case.toml")
    open(stream->TOML.print(stream,config),final_path,"w")
    result=(;source=file,output_directory=output,final_config_path=final_path,stopped=stopped[],history,runs,queries)
    open(stream->TOML.print(stream,_report_metadata(result)),joinpath(output,"batch.toml"),"w")
    result
end
