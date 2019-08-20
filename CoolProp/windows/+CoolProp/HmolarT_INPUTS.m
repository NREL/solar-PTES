function v = HmolarT_INPUTS()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 110);
  end
  v = vInitialized;
end
