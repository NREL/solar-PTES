function v = iDelta()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 26);
  end
  v = vInitialized;
end
