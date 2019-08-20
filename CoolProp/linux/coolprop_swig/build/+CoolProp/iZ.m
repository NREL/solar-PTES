function v = iZ()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 62);
  end
  v = vInitialized;
end
