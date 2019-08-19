function v = DmolarSmolar_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 129);
  end
  v = vInitialized;
end
