function v = iHelmholtzmolar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 35);
  end
  v = vInitialized;
end
