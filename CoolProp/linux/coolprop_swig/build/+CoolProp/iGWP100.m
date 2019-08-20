function v = iGWP100()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 68);
  end
  v = vInitialized;
end
