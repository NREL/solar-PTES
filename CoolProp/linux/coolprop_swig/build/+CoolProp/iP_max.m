function v = iP_max()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 19);
  end
  v = vInitialized;
end
