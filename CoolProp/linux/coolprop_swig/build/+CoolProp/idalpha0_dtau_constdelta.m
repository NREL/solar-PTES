function v = idalpha0_dtau_constdelta()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 56);
  end
  v = vInitialized;
end
