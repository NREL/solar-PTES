function v = HEOS_BACKEND_MIX()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 147);
  end
  v = vInitialized;
end
