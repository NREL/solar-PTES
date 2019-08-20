function v = idBvirial_dT()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 62);
  end
  v = vInitialized;
end
