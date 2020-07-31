function v = PUmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 121);
  end
  v = vInitialized;
end
