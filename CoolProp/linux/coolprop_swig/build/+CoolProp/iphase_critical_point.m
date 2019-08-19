function v = iphase_critical_point()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 80);
  end
  v = vInitialized;
end
