function v = iHelmholtzmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 45);
  end
  v = vInitialized;
end
