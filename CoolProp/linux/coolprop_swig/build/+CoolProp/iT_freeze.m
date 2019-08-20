function v = iT_freeze()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 66);
  end
  v = vInitialized;
end
