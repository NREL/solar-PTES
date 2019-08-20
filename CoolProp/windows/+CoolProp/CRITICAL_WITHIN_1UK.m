function v = CRITICAL_WITHIN_1UK()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 159);
  end
  v = vInitialized;
end
