function v = iGWP500()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 69);
  end
  v = vInitialized;
end
