function v = ialphar()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 52);
  end
  v = vInitialized;
end
