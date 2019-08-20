function v = iUmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 42);
  end
  v = vInitialized;
end
