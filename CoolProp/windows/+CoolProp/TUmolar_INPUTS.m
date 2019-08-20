function v = TUmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 114);
  end
  v = vInitialized;
end
