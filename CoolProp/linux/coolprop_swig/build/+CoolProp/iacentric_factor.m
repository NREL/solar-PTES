function v = iacentric_factor()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 6);
  end
  v = vInitialized;
end
