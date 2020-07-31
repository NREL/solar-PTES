function v = SmassUmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 124);
  end
  v = vInitialized;
end
