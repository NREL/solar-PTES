function v = TTSE_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 140);
  end
  v = vInitialized;
end
