function value = get_option_field(opts, field_name, default_value)
%GET_OPTION_FIELD Read an option struct field with a default fallback.

if isfield(opts, field_name)
    value = opts.(field_name);
else
    value = default_value;
end
end
