function v = TREND_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 139);
  end
  v = vInitialized;
end
