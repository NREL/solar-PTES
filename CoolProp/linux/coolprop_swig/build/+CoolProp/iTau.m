function v = iTau()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 25);
  end
  v = vInitialized;
end
