function v = idCvirial_dT()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 61);
  end
  v = vInitialized;
end
