function v = iT_min()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 17);
  end
  v = vInitialized;
end
