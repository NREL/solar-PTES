function v = SmolarUmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 125);
  end
  v = vInitialized;
end
