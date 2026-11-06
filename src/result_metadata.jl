# TOML has no null or Symbol/NamedTuple type. Keep scientific results in
# ordinary tables, explicitly spelling out an unavailable value.
function _result_metadata(value)
    if value isa NamedTuple || value isa AbstractDict
        # A missing optional reading ID must remain absent on settings reload;
        # the literal string "not_available" can itself be a valid input ID.
        return Dict(string(key)=>_result_metadata(item) for (key,item) in pairs(value)
            if !(string(key)=="observation_id" && item===nothing))
    elseif value isa Tuple || value isa AbstractVector
        return [_result_metadata(item) for item in value]
    elseif value isa Symbol
        return String(value)
    elseif value === nothing
        return "not_available"
    end
    return value
end

