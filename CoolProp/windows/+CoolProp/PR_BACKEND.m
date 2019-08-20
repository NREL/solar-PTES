function v = PR_BACKEND()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 156);
  end
  v = vInitialized;
end
