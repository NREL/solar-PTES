function v = INCOMP_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 137);
  end
  v = vInitialized;
end
