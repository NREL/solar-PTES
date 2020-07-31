function v = PR_BACKEND_FAMILY()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 143);
  end
  v = vInitialized;
end
