function v = CONFIGURATION_NOT_DEFINED_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 174);
  end
  v = vInitialized;
end
