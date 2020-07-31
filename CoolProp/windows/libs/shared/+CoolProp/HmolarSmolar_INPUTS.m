function v = HmolarSmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 125);
  end
  v = vInitialized;
end
