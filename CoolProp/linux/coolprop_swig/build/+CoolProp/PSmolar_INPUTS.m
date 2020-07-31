function v = PSmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 119);
  end
  v = vInitialized;
end
