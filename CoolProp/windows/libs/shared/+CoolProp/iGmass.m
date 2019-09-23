function v = iGmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 44);
  end
  v = vInitialized;
end
