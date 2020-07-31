function v = CRITICAL_SPLINES_ENABLED()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 134);
  end
  v = vInitialized;
end
