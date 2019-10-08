function v = iGWP100()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 70);
  end
  v = vInitialized;
end
