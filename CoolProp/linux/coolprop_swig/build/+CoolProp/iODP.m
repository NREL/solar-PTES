function v = iODP()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 73);
  end
  v = vInitialized;
end
