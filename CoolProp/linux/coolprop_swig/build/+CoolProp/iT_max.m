function v = iT_max()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 18);
  end
  v = vInitialized;
end
