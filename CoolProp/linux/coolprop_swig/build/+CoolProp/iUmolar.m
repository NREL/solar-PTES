function v = iUmolar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 33);
  end
  v = vInitialized;
end
