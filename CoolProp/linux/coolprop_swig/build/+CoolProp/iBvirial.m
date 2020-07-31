function v = iBvirial()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 58);
  end
  v = vInitialized;
end
