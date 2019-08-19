function v = iphase_supercritical_gas()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 78);
  end
  v = vInitialized;
end
