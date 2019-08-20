function v = INVALID_PARAMETER()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 3);
  end
  v = vInitialized;
end
