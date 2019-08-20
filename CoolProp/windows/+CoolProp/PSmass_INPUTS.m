function v = PSmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 120);
  end
  v = vInitialized;
end
