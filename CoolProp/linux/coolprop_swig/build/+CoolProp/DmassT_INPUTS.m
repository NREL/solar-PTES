function v = DmassT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 106);
  end
  v = vInitialized;
end
