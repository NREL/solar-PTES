function v = CONFIGURATION_ENDOFLIST_TYPE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 179);
  end
  v = vInitialized;
end
