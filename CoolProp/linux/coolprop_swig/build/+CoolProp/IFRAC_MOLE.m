function v = IFRAC_MOLE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 86);
  end
  v = vInitialized;
end
