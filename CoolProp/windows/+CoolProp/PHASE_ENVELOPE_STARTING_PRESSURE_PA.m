function v = PHASE_ENVELOPE_STARTING_PRESSURE_PA()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = CoolPropMATLAB_wrap(0, 171);
  end
  v = vInitialized;
end
