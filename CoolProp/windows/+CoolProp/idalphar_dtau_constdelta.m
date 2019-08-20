function v = idalphar_dtau_constdelta()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 55);
  end
  v = vInitialized;
end
