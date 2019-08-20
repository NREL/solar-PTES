function v = iundefined_parameter()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 75);
  end
  v = vInitialized;
end
