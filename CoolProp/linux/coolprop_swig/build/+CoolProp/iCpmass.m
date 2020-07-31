function v = iCpmass()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 39);
  end
  v = vInitialized;
end
