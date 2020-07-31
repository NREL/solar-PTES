function v = QSmass_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 100);
  end
  v = vInitialized;
end
