function v = iDmolar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 27);
  end
  v = vInitialized;
end
