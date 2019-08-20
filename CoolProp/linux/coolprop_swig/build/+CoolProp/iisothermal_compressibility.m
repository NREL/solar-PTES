function v = iisothermal_compressibility()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 49);
  end
  v = vInitialized;
end
