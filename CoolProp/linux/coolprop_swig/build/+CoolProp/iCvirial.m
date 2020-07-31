function v = iCvirial()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 59);
  end
  v = vInitialized;
end
