function v = iphase_supercritical_gas()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 80);
  end
  v = vInitialized;
end
