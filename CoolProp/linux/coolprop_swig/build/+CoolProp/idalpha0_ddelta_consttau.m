function v = idalpha0_ddelta_consttau()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 57);
  end
  v = vInitialized;
end
