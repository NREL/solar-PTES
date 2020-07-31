function v = iHH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 71);
  end
  v = vInitialized;
end
