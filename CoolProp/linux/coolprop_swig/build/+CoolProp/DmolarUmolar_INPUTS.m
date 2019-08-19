function v = DmolarUmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 131);
  end
  v = vInitialized;
end
