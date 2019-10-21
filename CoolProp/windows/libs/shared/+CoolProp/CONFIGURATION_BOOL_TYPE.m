function v = CONFIGURATION_BOOL_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 175);
  end
  v = vInitialized;
end
