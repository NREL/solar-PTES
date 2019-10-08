function v = TREND_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 152);
  end
  v = vInitialized;
end
