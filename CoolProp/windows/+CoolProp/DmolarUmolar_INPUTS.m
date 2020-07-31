function v = DmolarUmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 133);
  end
  v = vInitialized;
end
