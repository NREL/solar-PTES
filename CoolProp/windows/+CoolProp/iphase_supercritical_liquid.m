function v = iphase_supercritical_liquid()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 81);
  end
  v = vInitialized;
end
