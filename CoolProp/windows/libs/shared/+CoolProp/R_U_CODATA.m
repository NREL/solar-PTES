function v = R_U_CODATA()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 172);
  end
  v = vInitialized;
end
