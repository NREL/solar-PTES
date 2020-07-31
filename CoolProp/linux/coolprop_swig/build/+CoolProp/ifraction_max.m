function v = ifraction_max()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 65);
  end
  v = vInitialized;
end
