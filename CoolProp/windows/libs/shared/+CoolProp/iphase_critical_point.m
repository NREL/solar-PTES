function v = iphase_critical_point()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 82);
  end
  v = vInitialized;
end
