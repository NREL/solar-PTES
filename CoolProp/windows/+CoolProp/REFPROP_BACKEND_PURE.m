function v = REFPROP_BACKEND_PURE()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 148);
  end
  v = vInitialized;
end
