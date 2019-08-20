function v = irhomass_critical()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 12);
  end
  v = vInitialized;
end
