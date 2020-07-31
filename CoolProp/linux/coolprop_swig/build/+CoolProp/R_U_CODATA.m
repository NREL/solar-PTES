function v = R_U_CODATA()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolProp_wrap(0, 145);
  end
  v = vInitialized;
end
