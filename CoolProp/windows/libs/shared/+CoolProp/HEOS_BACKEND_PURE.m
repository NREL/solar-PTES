function v = HEOS_BACKEND_PURE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 146);
  end
  v = vInitialized;
end
