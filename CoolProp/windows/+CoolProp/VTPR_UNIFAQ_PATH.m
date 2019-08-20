function v = VTPR_UNIFAQ_PATH()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 173);
  end
  v = vInitialized;
end
