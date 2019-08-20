function v = iphase_supercritical_liquid()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 79);
  end
  v = vInitialized;
end
