const _CLI_HELP = """
MaoyeLoc — 用地震波到时及到时差估计震源位置和发震时刻

用法：
  myloc CASE.toml
  myloc --verbose CASE.toml
  myloc --quiet CASE.toml
  myloc --diagnostics CASE.toml
  myloc --check CASE.toml
  myloc --show-config CASE.toml
  myloc --template NEW_CASE.toml
  myloc --batch BATCH.toml
  myloc --version
  myloc --help

选项：
  --check        读取模型与输入并做静态检查，不定位、不写出或清理结果。
  --show-config  将实际生效参数及显式填写的参数名写到标准输出，不运行定位。
  --template     写入带注释的待填写 TOML 模板；拒绝覆盖已有路径。
  --verbose      显示逐事件/双差连通分量进度与最终详单。
  --quiet        只打印最终阶段状态、观测计数与输出位置，不显示进度。
  --diagnostics  输出观测残差、权重、选择及迭代等诊断文件。
  --batch        在同一进程中按 TOML steps 运行批处理。
  --version      显示程序版本。
  --help, -h     显示本帮助。

CASE 中 task 可为 absolute、relative 或 absolute_then_relative。
分别表示绝对定位、双差重定位，或先绝对定位再双差重定位。
CT 为目录拾取形成的到时差，CC 为外部波形互相关测得的到时差。
双差连通分量是通过差分观测直接或间接相连的一组事件。
按 [files]、[model]、[solver]、[relative]、[uncertainty] 组织输入与参数。
绝对定位只有一条路径：经验加权、SVD 改正和步长控制。
相对定位只有一条路径：分阶段改权、列缩放和阻尼。详见用户手册。
文件路径相对于 CASE 所在目录；内部单位为 km、km/s、s，x东/y北/z下。
必须明确 UTC 时间参考、坐标及垂直基准；模板不提供通用速度模型。
静态检查通过不代表全部射线可达或全部事件收敛，结果须结合实际 status。
默认按真实完成量报告有限次进度；阶段完成不等于数值收敛。
可选 [depth_reporting] 的近站选站规则及报告深度含义见用户手册。
参数细节见 docs/用户手册.md，字段定义见 docs/数据约定.md。
"""

const _CLI_CASE_TEMPLATE = raw"""
# MaoyeLoc 原生配置模板。所有空字符串与空模型数组均须先填写。
# 本模板没有通用速度模型；模型类型须按实际物理问题选择。
# 路径相对于本 TOML 所在目录；事件/台站 ID 按字符串原样保留。
task = "absolute_then_relative" # absolute / relative / absolute_then_relative
time_reference_utc = ""       # 必填：共同 UTC 参考，ISO8601 且以 Z 结尾
coordinates = ""             # 必填：坐标定义；内部 x东/y北/z下，长度 km
vertical_reference = ""      # 必填：源、站和模型层界共用的垂直基准
output_directory = "results" # 本次结果目录；复用目录会更新所请求阶段的产物
# reference_event_ids = ["EVENT_ID"] # 固定参考事件的位置与起时；此键位于顶层

[files]
stations = ""                 # station_id,x_km,y_km,z_km
events = ""                   # event_id,x_km,y_km,z_km,origin_s
picks = ""                    # event_id,station_id,phase,arrival_s,sigma_s,quality_code
# 绝对定位须提供quality_code，或显式配置absolute_weighting.default_quality。
differentials = ""            # event_i,event_j,station_id,phase,dt_s,sigma_s,source
# arrival_s、origin_s 相对上述共同参考；dt_s = arrival_i - arrival_j。
# 若仅做 absolute，删除 differentials；relative 的观测来源须显式给出。

[model]
# 这里展示分层常速模型的字段结构；其他已有模型类型见用户手册。
kind = "constant_velocity_layers"
layer_top_km = []              # 必填：严格递增的层顶，km，与源站同一基准
vp_km_s = []                   # 必填：各层 P 速度，km/s
vs_km_s = []                   # 必填：各层 S 速度，km/s；每层 Vp > Vs > 0

[solver]
method = "hypoinverse"        # 唯一绝对定位路径：经验加权、SVD 改正、步长控制
# 具体生效权、谱控制与停止参数由 --show-config 列出。
# depth_min_km = ...           # 默认模型顶面，km；需要时显式声明
# depth_max_km = ...           # 默认无上限，km
# max_iterations = 50
# location_tolerance_km = 1e-6
# time_tolerance_s = 1e-7

[relative]
# 唯一相对定位路径：分阶段按类别、距离和残差改权，列缩放后加阻尼。
# 小系统用同一阻尼问题的密集 SVD，大系统用 LSQR。不自动固定参考事件。
# use_absolute_picks = false   # true 表示同时使用绝对到时；检查共享误差假设
# initialization = "catalog" # 保留目录初值；不自动改成质心
# damp = 5.0                 # 覆盖默认各阶段阻尼；作用在列缩放后的改正量上

[uncertainty]
# 仅在误差假设已核实时显式请求；不默认宣称观测独立。
# assume_independent_picks = false
# assume_independent_relative_observations = false
# relative_pairs = [["EVENT_I", "EVENT_J"]]

[execution]
tasks = 1                    # 独立绝对定位事件/完整双差连通分量的并发上限
replay_tasks = 1             # 沿同一观测误差方向扰动后完整重算的并发上限，最多4
# tasks与replay_tasks不可同时大于1；可执行程序用 --julia-args --threads=N 设置线程。

[output]
diagnostics = false          # 默认保存定位目录及已请求误差，true增加诊断文件
detail = "full"              # 诊断开启后：full包含历史，core省略逐次历史和谱方向

# 可选报告深度；以下整段默认不启用，须先核实源台及模型的垂直基准。
# [depth_reporting]
# model_datum = ""            # 必填：nearest5用所选站平均高程定义基准面；geoid用所声明的大地水准面
# coordinate_datum = ""       # 必填：计算z采用模型坐标(model)还是以声明的大地水准面为零面(geoid)
# report_datum = ""           # 必填：报告深度相对于model基准面还是所声明的大地水准面geoid
# station_selection = "nearest5" # 默认不设距离/高程门槛；可显式改为 nearby5
# nearby5：距离≤100 km且高程≥−0.120 km，最多五个不同站点；无合格近站
# 时回退到全部候选中最近的高程已知记录，允许超距/低高程，未知不补零。
# model坐标按所选站平均高程换算深度时，须另给station_elevations_km高程表(km)。
# 此处不计算大地水准面起伏，也不将椭球高自动换为物理高程。
# 不移动计算源台；报告深度及使用该深度的Hypo71质量码可能变化。
# nearby5另导出datum_stations.csv记录逐站取舍；完整映射写法见手册。
"""

function _cli_status_counts(rows)
    counts = Dict{String,Int}()
    for row in rows
        status = String(row.status)
        counts[status] = get(counts,status,0)+1
    end
    isempty(counts) ? "(无)" : join(("$status=$number" for (status,number) in sort!(collect(counts);by=first)),", ")
end

function _cli_observation_summary(label,rows)
    used = count(row->row.used,rows)
    println(label, ": 输入 ",length(rows),"，使用 ",used,"，未用 ",length(rows)-used)
end

function _cli_print_summary(result;verbose=false)
    if result.absolute !== nothing
        absolute = result.absolute
        println("绝对定位事件 ",length(absolute.events),"：",_cli_status_counts(absolute.events))
        _cli_observation_summary("绝对到时",absolute.residuals)
        println("绝对目录筛除：",count(event->!get(event,:catalog_eligible,true),absolute.events))
        if verbose
            for event in absolute.events
                println("  ",event.event_id," absolute: ",event.status,"; used ",event.used_count,"/",event.input_count,
                    get(event,:catalog_eligible,true) ? "" : "; catalog excluded: "*String(event.catalog_exclusion_reason))
            end
        end
    end
    if result.relative !== nothing
        relative = result.relative
        println("相对定位事件 ",length(relative.events),"：",_cli_status_counts(relative.events))
        println("双差连通分量 ",length(relative.components),"：",_cli_status_counts(relative.components))
        _cli_observation_summary("差分观测",relative.residuals)
        _cli_observation_summary("相对阶段绝对到时",relative.absolute_residuals)
        not_entered = count(event->event.status==:not_entered_relative,relative.events)
        selected_out = count(event->get(event,:selection_status,:selected)!=:selected,relative.events)
        exited = count(event->event.status in (:airquake_removed,:disconnected_after_airquake),relative.events)
        println("未进入相对阶段（绝对定位未通过）：",not_entered,
            "；选择筛除：",selected_out,"；试步深度越界（airquake）相关退出：",exited)
        if verbose
            for group in relative.components
                println("  双差连通分量 ",group.component_id,": ",group.status,"; 整体位置与共同起时的约束（datum） ",group.datum,
                    isempty(group.reference_event_ids) ? "" : "; reference "*join(group.reference_event_ids,","),
                    "; used ",group.used_count,"/",group.input_count)
            end
            for event in relative.events
                println("  ",event.event_id," relative: ",event.status,"; used ",event.used_count,"/",event.input_count)
            end
        end
    end
    println("结果目录：",result.output_directory)
    println("各 status 分别保留收敛、阶段完成、固定、失败和未用含义；详见 events.csv 与本次运行的 run_summary.toml。")
end

function _cli_progress_reporter(verbose)
    last_reported = Dict{Symbol,Int}()
    return function(info)
        stage = Symbol(info.stage)
        if get(info,:kind,nothing)==:phase
            elapsed=get(info,:elapsed_seconds,nothing)
            println(stderr,"[",stage,"] ",info.status,
                elapsed===nothing ? "" : " ("*string(round(elapsed;digits=3))*" s)")
            return nothing
        end
        done,total = info.processed,info.total
        replay=stage in (:relative_replay_baseline,:relative_replay_directions)
        stride = max(replay ? 1 : 100,cld(total,10))
        if verbose || done==1 || done==total || get(info,:status,nothing) in (:started,:stopped,:completed) || done-get(last_reported,stage,0)>=stride
            event_id,component_id = get(info,:event_id,nothing),get(info,:component_id,nothing)
            identity = event_id!==nothing ? " event "*String(event_id) :
                component_id!==nothing ? " component "*string(component_id) : ""
            println(stderr,"[",stage,"] ",done,"/",total,verbose ? identity*"; "*String(info.status) : "")
            last_reported[stage] = done
        end
        nothing
    end
end

function cli_main(args=ARGS)
    if isempty(args) || args in (["--help"],["-h"])
        print(_CLI_HELP)
        return
    end
    if args==["--version"]
        println("MaoyeLoc ", _APP_VERSION)
        return
    elseif first(args)=="--batch"
        length(args)==2 || throw(ArgumentError("用法：myloc --batch BATCH.toml"))
        result=run_batch(last(args))
        println("批处理结果目录：",result.output_directory)
        return
    end
    mode = startswith(first(args),"-") ? first(args) : "run"
    mode in ("run","--check","--show-config","--template","--verbose","--quiet","--diagnostics") ||
        throw(ArgumentError("未知选项 $(first(args))；使用 --help 查看用法。"))
    length(args)==(mode=="run" ? 1 : 2) ||
        throw(ArgumentError("此命令需要一个配置路径；使用 --help 查看用法。"))
    path = last(args)
    if mode=="--template"
        (ispath(path) || islink(path)) && throw(ArgumentError("模板目标已存在，拒绝覆盖：$path"))
        write(path,_CLI_CASE_TEMPLATE)
        println("已写入待填写模板：",abspath(path))
        println("请先填写时空定义、文件路径和模型，再运行 --check。")
        return
    end
    prepared = if mode=="--diagnostics"
        config=TOML.parsefile(path)
        config["output"]=merge(get(config,"output",Dict()),Dict("diagnostics"=>true))
        _prepare_case(config;base_directory=dirname(abspath(path)),config_path=abspath(path))
    else
        prepare_case(path)
    end
    if mode=="--check"
        println("# 静态输入检查；未执行定位，不保证每条射线可达或每个事件收敛。")
        TOML.print(stdout,check_case(prepared))
        return
    elseif mode=="--show-config"
        TOML.print(stdout,effective_config(prepared))
        return
    end
    quiet = mode=="--quiet"
    verbose = mode=="--verbose"
    if !quiet
        println("任务：",prepared.task,"；事件 ",length(prepared.events),"；台站 ",length(prepared.stations),
            "；绝对到时 ",prepared.picks===nothing ? "未读入" : length(prepared.picks),
            "；差分 ",prepared.differences===nothing ? "未读入/按配置生成" : length(prepared.differences),
            "；请求并发任务 ",prepared.tasks,"（当前线程最多 ",min(prepared.tasks,Threads.nthreads()),"）",
            "；请求观测扰动后完整重算并发 ",prepared.replay_tasks,"（最多 ",min(prepared.replay_tasks,4,Threads.nthreads()),"）")
    end
    result = run_case(prepared;progress=quiet ? nothing : _cli_progress_reporter(verbose))
    _cli_print_summary(result;verbose)
end

# PackageCompiler and the source scripts enter the same command path.
function julia_main()::Cint
    try
        cli_main(ARGS)
        return 0
    catch exception
        exception isa InterruptException && rethrow()
        print(stderr,"MaoyeLoc: ")
        showerror(stderr,exception)
        println(stderr)
        return 1
    end
end
