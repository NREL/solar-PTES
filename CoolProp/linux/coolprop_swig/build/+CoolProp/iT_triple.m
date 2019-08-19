function v = iT_triple()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 15);
  end
  v = vInitialized;
end
