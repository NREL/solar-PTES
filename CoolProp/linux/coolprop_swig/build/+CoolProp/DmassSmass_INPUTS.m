function v = DmassSmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 128);
  end
  v = vInitialized;
end
