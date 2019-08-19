function v = DONT_CHECK_PROPERTY_LIMITS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 142);
  end
  v = vInitialized;
end
