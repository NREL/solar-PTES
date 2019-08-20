function v = iphase_unknown()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 83);
  end
  v = vInitialized;
end
