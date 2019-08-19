function v = CONFIGURATION_INTEGER_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 149);
  end
  v = vInitialized;
end
