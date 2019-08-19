function v = iphase_liquid()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 76);
  end
  v = vInitialized;
end
