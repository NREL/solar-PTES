function v = iSmolar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 29);
  end
  v = vInitialized;
end
