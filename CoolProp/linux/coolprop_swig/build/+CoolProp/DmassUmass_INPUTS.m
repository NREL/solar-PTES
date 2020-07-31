function v = DmassUmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 130);
  end
  v = vInitialized;
end
