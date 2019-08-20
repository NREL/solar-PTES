function v = idalphar_dtau_constdelta()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 53);
  end
  v = vInitialized;
end
