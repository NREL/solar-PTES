function v = iGmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 43);
  end
  v = vInitialized;
end
