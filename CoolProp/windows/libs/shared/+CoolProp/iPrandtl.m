function v = iPrandtl()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 49);
  end
  v = vInitialized;
end
