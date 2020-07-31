function v = INVALID_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 134);
  end
  v = vInitialized;
end
