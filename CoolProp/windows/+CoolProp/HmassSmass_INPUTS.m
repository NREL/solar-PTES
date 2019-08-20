function v = HmassSmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 124);
  end
  v = vInitialized;
end
