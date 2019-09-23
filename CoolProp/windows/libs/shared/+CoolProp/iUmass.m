function v = iUmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 43);
  end
  v = vInitialized;
end
