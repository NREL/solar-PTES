function v = REFPROP_USE_GERG()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 140);
  end
  v = vInitialized;
end
