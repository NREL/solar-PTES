function v = IFRAC_MASS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 85);
  end
  v = vInitialized;
end
