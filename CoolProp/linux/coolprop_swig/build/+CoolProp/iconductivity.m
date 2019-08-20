function v = iconductivity()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 45);
  end
  v = vInitialized;
end
