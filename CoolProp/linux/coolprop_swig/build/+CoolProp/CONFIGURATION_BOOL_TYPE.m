function v = CONFIGURATION_BOOL_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 147);
  end
  v = vInitialized;
end
