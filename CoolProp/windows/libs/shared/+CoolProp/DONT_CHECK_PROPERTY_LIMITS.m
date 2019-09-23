function v = DONT_CHECK_PROPERTY_LIMITS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 169);
  end
  v = vInitialized;
end
