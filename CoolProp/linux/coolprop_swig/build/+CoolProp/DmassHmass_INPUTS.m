function v = DmassHmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 126);
  end
  v = vInitialized;
end
