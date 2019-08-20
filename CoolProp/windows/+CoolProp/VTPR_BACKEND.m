function v = VTPR_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 157);
  end
  v = vInitialized;
end
