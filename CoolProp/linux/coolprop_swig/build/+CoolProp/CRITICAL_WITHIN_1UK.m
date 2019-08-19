function v = CRITICAL_WITHIN_1UK()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 133);
  end
  v = vInitialized;
end
