function v = CONFIGURATION_DOUBLE_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 176);
  end
  v = vInitialized;
end
