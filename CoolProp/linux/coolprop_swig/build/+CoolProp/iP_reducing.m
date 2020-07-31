function v = iP_reducing()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 14);
  end
  v = vInitialized;
end
