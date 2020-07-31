function v = iphase_gas()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 81);
  end
  v = vInitialized;
end
