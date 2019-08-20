function v = IFRAC_MOLE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 88);
  end
  v = vInitialized;
end
