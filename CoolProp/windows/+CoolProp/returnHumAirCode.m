function varargout = returnHumAirCode(varargin)
  [varargout{1:max(1,nargout)}] = CoolPropMATLAB_wrap(461,varargin{:});
end
