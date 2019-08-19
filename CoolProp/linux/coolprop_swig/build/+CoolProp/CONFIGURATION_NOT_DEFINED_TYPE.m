function v = CONFIGURATION_NOT_DEFINED_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 146);
  end
  v = vInitialized;
end
