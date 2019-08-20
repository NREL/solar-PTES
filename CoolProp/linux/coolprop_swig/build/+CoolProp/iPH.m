function v = iPH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 72);
  end
  v = vInitialized;
end
