function v = ialpha0()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 55);
  end
  v = vInitialized;
end
