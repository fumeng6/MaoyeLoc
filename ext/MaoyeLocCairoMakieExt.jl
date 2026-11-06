module MaoyeLocCairoMakieExt

using MaoyeLoc
using CairoMakie
using TOML

function plot_rows(directory,name,required)
    rows=Dict{String,String}[]
    MaoyeLoc._foreach_csv_row(joinpath(directory,name),required;
            allow_extra=true,allow_empty=true,strip_unquoted=false) do line,fields,header
        push!(rows,Dict(zip(header,fields)))
    end
    rows
end

numbers(rows,key)=Float64[parse(Float64,row[key]) for row in rows]

function MaoyeLoc._render_location_plot(directory::AbstractString)
    CairoMakie.activate!()
    config=TOML.parsefile(joinpath(directory,"plot.toml"))
    events=plot_rows(directory,"relocated.csv",("event_id","x_km","y_km","depth_km","origin_s","status"))
    stations=unique(row->row["station_id"],plot_rows(directory,"stations.csv",("station_id","x_km","y_km")))
    frame=config["frame"]
    geographic=frame["kind"]=="wgs84_local"
    xkey,ykey=geographic ? ("longitude_deg","latitude_deg") : ("x_km","y_km")
    xlabel,ylabel=geographic ? ("Longitude (°)","Latitude (°)") : ("East (km)","North (km)")
    points=vcat(events,stations)
    if isempty(points)
        throw(ArgumentError("No event or station coordinates available for plotting."))
    end
    x,y=numbers(points,xkey),numbers(points,ykey)
    cx,cy=(minimum(x)+maximum(x))/2,(minimum(y)+maximum(y))/2
    xscale=geographic ? max(cosd(frame["latitude_deg"]),1e-6) : 1.
    span=1.12max((maximum(x)-minimum(x))*xscale,maximum(y)-minimum(y),geographic ? .001 : .1)
    xlimits=(cx-span/(2xscale),cx+span/(2xscale))
    ylimits=(cy-span/2,cy+span/2)
    depth=numbers(events,"depth_km")
    if isempty(depth)
        colorlimits=(0.,1.)
        depthlimits=(0.,1.)
    else
        lower,upper=extrema(depth)
        padding=max(.05max(upper-lower,abs(upper)),.01)
        colorlimits=lower==upper ? (lower-padding,upper+padding) : (lower,upper)
        depthlimits=(min(0.,lower-padding),upper+padding)
    end
    colormap=Reverse(:inferno)
    statuses=sort!(unique(row["status"] for row in events))
    markers=Dict("converged"=>:circle,"reference_fixed"=>:diamond,"fixed_all"=>:diamond,
        "no_observations"=>:xcross,"failed"=>:xcross,"supplied_catalog"=>:circle)
    fallback=(:rect,:cross,:utriangle,:pentagon,:star5,:dtriangle)
    for (index,status) in enumerate(statuses)
        get!(markers,status,fallback[mod1(index,length(fallback))])
    end
    theme=Theme(fontsize=12,textcolor="#27343D",
        Axis=(xlabelsize=12,ylabelsize=12,titlesize=12,titlealign=:left,
            xticklabelsize=11,yticklabelsize=11,backgroundcolor="#F4F5F6",
            xgridcolor=(:black,.08),ygridcolor=(:black,.08),
            spinewidth=.7,xtickwidth=.7,ytickwidth=.7,xticksize=4,yticksize=4,
            xticks=WilkinsonTicks(5;k_min=4,k_max=6),yticks=WilkinsonTicks(5;k_min=4,k_max=6)),
        Legend=(framevisible=false,labelsize=11,patchsize=(14,14),rowgap=4),
        Colorbar=(labelsize=12,ticklabelsize=11,spinewidth=.5))
    with_theme(theme) do
        figure=Figure(figure_padding=20)
        Label(figure[0,1:2],"Event locations and cross-sections · N = $(length(events))";
            fontsize=18,font=:bold,halign=:left,tellwidth=false)
        horizontal=Axis(figure[1,1];width=480,height=480,xlabel,ylabel,
            title="(a) Horizontal distribution",aspect=AxisAspect(1))
        north_depth=Axis(figure[1,2];width=240,height=480,xlabel="Depth (km)",
            title="(b) Depth–$(geographic ? "latitude" : "north")",yticklabelsvisible=false)
        east_depth=Axis(figure[2,1];width=480,height=230,xlabel,ylabel="Depth (km)",
            title="(c) $(geographic ? "Longitude" : "East")–depth",valign=:top)
        marker_size=clamp(22/sqrt(max(length(events),1)),2.,7.)
        for status in statuses
            selected=[row for row in events if row["status"]==status]
            z=numbers(selected,"depth_km")
            for (axis,xvalues,yvalues) in ((horizontal,numbers(selected,xkey),numbers(selected,ykey)),
                    (north_depth,z,numbers(selected,ykey)),(east_depth,numbers(selected,xkey),z))
                scatter!(axis,xvalues,yvalues;color=z,colormap,colorrange=colorlimits,
                    marker=markers[status],markersize=marker_size,
                    strokecolor=("#27343D",.55),strokewidth=length(events)<=200 ? .35 : 0.)
            end
        end
        if !isempty(stations)
            scatter!(horizontal,numbers(stations,xkey),numbers(stations,ykey);
                marker=:utriangle,markersize=9,color=:transparent,strokecolor="#45545E",strokewidth=1)
        end
        xlims!(horizontal,xlimits...); ylims!(horizontal,ylimits...)
        xlims!(north_depth,depthlimits...); ylims!(north_depth,ylimits...)
        xlims!(east_depth,xlimits...); ylims!(east_depth,reverse(depthlimits)...)
        linkyaxes!(horizontal,north_depth)
        linkxaxes!(horizontal,east_depth)
        if isempty(events)
            for axis in (horizontal,north_depth,east_depth)
                text!(axis,.5,.5;text="No result events",space=:relative,align=(:center,:center),color="#65737E")
            end
        end

        times=sort!(numbers([row for row in events if !isempty(row["origin_s"])],"origin_s"))
        duration=isempty(times) ? 0. : last(times)-first(times)
        unit,label=duration>=86400 ? (86400.,"d") : duration>=3600 ? (3600.,"h") :
            duration>=60 ? (60.,"min") : (1.,"s")
        timeline=Axis(figure[2,2];width=240,height=230,valign=:top,
            title="(d) Count · n=$(length(times))/$(length(events))",xlabel="Elapsed time ($label)",ylabel="Events",
            xticks=WilkinsonTicks(4;k_min=3,k_max=5),
            yticks=length(times)<=10 ? collect(0:max(1,cld(length(times),5)):max(length(times),1)) :
                WilkinsonTicks(5;k_min=4,k_max=6),backgroundcolor=:white)
        if isempty(times)
            text!(timeline,.5,.5;text="No origin times",space=:relative,align=(:center,:center),color="#65737E")
            xlims!(timeline,0,1); ylims!(timeline,0,1)
        else
            elapsed=(times.-first(times))./unit
            stairs!(timeline,[0.;elapsed],[0;collect(eachindex(times))];
                step=:post,color="#416B85",linewidth=1.6)
            xlims!(timeline,0,max(1.04duration/unit,.1)); ylims!(timeline,0,1.05length(times))
        end
        if !isempty(events)
            scale=figure[3,2]=GridLayout()
            Label(scale[1,1],"Depth (km)";fontsize=12,halign=:right)
            Colorbar(scale[1,2];colormap,limits=colorlimits,
                vertical=false,flipaxis=false,width=160,height=10,tellwidth=false,
                ticks=WilkinsonTicks(4;k_min=3,k_max=5))
            colgap!(scale,8)
        end

        elements=[MarkerElement(marker=markers[status],color="#45545E",markersize=6) for status in statuses]
        labels=["$status ($(count(row->row["status"]==status,events)))" for status in statuses]
        if !isempty(stations)
            push!(elements,MarkerElement(marker=:utriangle,color=:transparent,strokecolor="#45545E",strokewidth=1))
            push!(labels,"Station ($(length(stations)))")
        end
        if !isempty(elements)
            Legend(figure[3,1],elements,labels;orientation=:horizontal,nbanks=max(1,cld(length(labels),4)),
                halign=:left,tellwidth=false,tellheight=true)
        end
        colgap!(figure.layout,18); rowgap!(figure.layout,14)
        resize_to_layout!(figure)
        save(joinpath(directory,"locations.png"),figure;px_per_unit=2)
        save(joinpath(directory,"locations.pdf"),figure)
    end
    directory
end

end
