function v = iisobaric_expansion_coefficient()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 52);
  end
  v = vInitialized;
end
