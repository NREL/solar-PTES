function v = CONFIGURATION_DOUBLE_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 148);
  end
  v = vInitialized;
end
